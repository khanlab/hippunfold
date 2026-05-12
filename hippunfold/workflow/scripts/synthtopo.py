#!/usr/bin/env python3
"""
Fully isolated SynthSeg inference script for hippocampus segmentation.
This script contains all necessary code and does not require local imports.

Author: Mahmoud Salman (mahmoud1yaser)
"""

import argparse
import os
import time
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import nibabel as nib
from typing import Sequence, Tuple, Optional
from types import GeneratorType as generator

from cornucopia import LoadTransform
from cornucopia.utils import warps  # original sampling ops


# utils.ensure_list equivalent
def _ensure_list(x, size=None, crop=True):
    if not isinstance(x, (list, tuple, range)):
        x = [x]
    elif not isinstance(x, list):
        x = list(x)
    if size and len(x) < size:
        x += x[-1:] * (size - len(x))
    if size and crop:
        x = x[:size]
    return x


# === Modules (subset of learn2synth.modules) ===
class ConvBlockBase(nn.Sequential):
    """Original ConvBlockBase (order auto-fix, norm channel inference)."""

    def __init__(
        self,
        ndim,
        in_channels,
        out_channels,
        opt_conv=None,
        activation="LeakyReLU",
        norm=None,
        dropout=False,
        order="cand",
    ):
        super().__init__()
        self.order = self._fix_order(order)
        conv_cls = getattr(nn, f"Conv{ndim}d")
        opt_conv = dict(opt_conv or {})
        opt_conv.setdefault("kernel_size", 3)
        opt_conv.setdefault("bias", True)
        opt_conv.setdefault("padding", "same")
        conv = conv_cls(in_channels, out_channels, **opt_conv)
        act = self._make_activation(activation)
        drop = self._make_dropout(dropout, ndim)
        norm_mod = self._make_norm(norm, ndim, conv, self.order)
        for ch in self.order:
            if ch == "n" and norm_mod is not None:
                self.add_module("norm", norm_mod)
            elif ch == "c":
                self.add_module("conv", conv)
            elif ch == "d" and drop is not None:
                self.add_module("dropout", drop)
            elif ch == "a" and act is not None:
                self.add_module("activation", act)

    @staticmethod
    def _fix_order(order):
        order = order.lower()
        for ch in "ncda":
            if ch not in order:
                order += ch
        return order

    @staticmethod
    def _make_activation(activation):
        if not activation:
            return None
        if isinstance(activation, str):
            activation = getattr(nn, activation)
        return activation() if isinstance(activation, type) else activation

    @staticmethod
    def _make_dropout(dropout, ndim):
        if not dropout:
            return None
        if isinstance(dropout, (int, float)):
            return getattr(nn, f"Dropout{ndim}d")(p=float(dropout))
        return dropout() if isinstance(dropout, type) else dropout

    @staticmethod
    def _make_norm(norm, ndim, conv, order):
        if not norm:
            return None
        if isinstance(norm, bool) and norm:
            norm = "batch"
        idx_n = order.index("n")
        idx_c = order.index("c")
        in_ch = conv.in_channels if idx_n < idx_c else conv.out_channels
        if isinstance(norm, str):
            if "instance" in norm.lower():
                norm_cls = getattr(nn, f"InstanceNorm{ndim}d")
                return norm_cls(in_ch)
            if "batch" in norm.lower():
                norm_cls = getattr(nn, f"BatchNorm{ndim}d")
                return norm_cls(in_ch)
            if "layer" in norm.lower():
                return nn.GroupNorm(in_ch, in_ch)
        return norm() if isinstance(norm, type) else norm


class ConvBlock(ConvBlockBase):
    def __init__(
        self,
        ndim,
        in_channels,
        out_channels,
        kernel_size=3,
        bias=True,
        activation="LeakyReLU",
        norm="instance",
        dropout=0,
        order="cand",
    ):
        super().__init__(
            ndim,
            in_channels,
            out_channels,
            opt_conv=dict(kernel_size=kernel_size, bias=bias, padding="same"),
            activation=activation,
            norm=norm,
            dropout=dropout,
            order=order,
        )


class ConvGroup(nn.Sequential):
    """Original ConvGroup logic (handles in/mid/out transitions, residual)."""

    def __init__(
        self,
        ndim,
        in_channels,
        mid_channels=None,
        out_channels=None,
        kernel_size=3,
        nb_conv=2,
        residual=False,
        activation="LeakyReLU",
        norm="instance",
        dropout=0,
        order="cand",
    ):
        self.residual = residual
        mid_channels = mid_channels or in_channels
        out_channels = out_channels or mid_channels
        nb_conv = (
            nb_conv - (in_channels != mid_channels) - (out_channels != mid_channels)
        )
        layers = []
        if in_channels != mid_channels:
            layers.append(
                ConvBlock(
                    ndim,
                    in_channels,
                    mid_channels,
                    kernel_size,
                    activation=activation,
                    norm=norm,
                    dropout=dropout,
                    order=order,
                )
            )
        for _ in range(nb_conv):
            layers.append(
                ConvBlock(
                    ndim,
                    mid_channels,
                    mid_channels,
                    kernel_size,
                    activation=activation,
                    norm=norm,
                    dropout=dropout,
                    order=order,
                )
            )
        if out_channels != mid_channels:
            layers.append(
                ConvBlock(
                    ndim,
                    mid_channels,
                    out_channels,
                    kernel_size,
                    activation=activation,
                    norm=norm,
                    dropout=dropout,
                    order=order,
                )
            )
        super().__init__(*layers)

    def forward(self, x):
        if self.residual:
            for layer in self:
                ident = x
                x = layer(x)
                if x.shape[1] == ident.shape[1]:
                    x = x + ident
            return x
        return super().forward(x)


class Downsample(nn.Module):
    """Original warps-based Downsample."""

    def __init__(self, factor=2, anchor="center"):
        super().__init__()
        self.factor = factor
        self.anchor = anchor

    def forward(self, x, shape=None):
        factor = None if shape else self.factor
        return warps.downsample(x, factor, shape, self.anchor)


class Upsample(nn.Module):
    """Original warps-based Upsample."""

    def __init__(self, factor=2, anchor="center"):
        super().__init__()
        self.factor = factor
        self.anchor = anchor

    def forward(self, x, shape=None):
        factor = None if shape else self.factor
        return warps.upsample(x, factor, shape, self.anchor)


class ConvBlockDown(nn.Sequential):
    """Original ConvBlockDown (downsample + conv group)."""

    def __init__(
        self,
        ndim,
        in_channels,
        out_channels,
        factor=2,
        kernel_size=3,
        activation="LeakyReLU",
        norm="instance",
        dropout=0,
        order="cand",
    ):
        super().__init__()
        self.downsample = Downsample(factor=factor)
        self.conv = ConvBlock(
            ndim,
            in_channels,
            out_channels,
            kernel_size=kernel_size,
            activation=activation,
            norm=norm,
            dropout=dropout,
            order=order,
        )

    def forward(self, x):
        return self.conv(self.downsample(x))


class ConvBlockUp(nn.Module):
    """Original ConvBlockUp (conv + upsample + skip combine)."""

    def __init__(
        self,
        ndim,
        in_channels,
        out_channels,
        factor=2,
        kernel_size=3,
        activation="LeakyReLU",
        norm="instance",
        dropout=0,
        order="cand",
        combine="cat",
    ):
        super().__init__()
        self.conv = ConvBlock(
            ndim,
            in_channels,
            out_channels,
            kernel_size=kernel_size,
            activation=activation,
            norm=norm,
            dropout=dropout,
            order=order,
        )
        self.upsample = Upsample(factor=factor)
        self.combine = combine

    def forward(self, x, skip=None):
        shape = None if skip is None else skip.shape[2:]
        x = self.conv(x)
        x = self.upsample(x, shape)
        if skip is not None:
            x = torch.cat([x, skip], dim=1) if self.combine == "cat" else x + skip
        return x


class EncoderBlock(nn.Sequential):
    def __init__(self, down, conv):
        super().__init__()
        self.down = down
        self.conv = conv

    def forward(self, x):
        return self.conv(self.down(x))


class DecoderBlock(nn.Sequential):
    def __init__(self, conv, up):
        super().__init__()
        self.conv = conv
        self.up = up

    def forward(self, x, skip=None):
        x = self.conv(x)
        return self.up(x, skip)


class UNet(nn.Module):
    """Original-style UNet (encoder/decoder pyramids with warps sampling)."""

    class defaults:
        nb_levels = 6
        nb_features = (16, 24, 32, 48, 64, 96, 128, 192, 256, 320)
        nb_conv = 2
        kernel_size = 3
        activation = "LeakyReLU"
        norm = "instance"
        dropout = 0
        residual = False
        factor = 2
        use_strides = False
        order = "cand"
        combine = "cat"

    def __init__(self, ndim, **kwargs):
        super().__init__()
        # load defaults
        for k, v in UNet.defaults.__dict__.items():
            if not k.startswith("_"):
                kwargs.setdefault(k, v)
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.ndim = ndim
        self.nb_features = _ensure_list(self.nb_features, self.nb_levels)
        self.in_channels = self.out_channels = self.nb_features[0]
        # encoder
        i, o = self.nb_features[0], self.nb_features[0]
        encoder = [self._conv_block(i, o)]
        for n in range(1, len(self.nb_features) - 1):
            i, o = self.nb_features[n - 1], self.nb_features[n]
            encoder.append(EncoderBlock(self._down_block(i, o), self._conv_block(o)))
        if self.nb_levels > 1:
            i, o = self.nb_features[-2], self.nb_features[-1]
            encoder.append(self._down_block(i, o))
        self.encoder = nn.Sequential(*encoder)
        # decoder
        decoder = []
        for n in range(len(self.nb_features) - 1):
            i, o = self.nb_features[-n - 1], self.nb_features[-n - 2]
            m = i
            if self.combine == "cat" and n > 0:
                i *= 2
            decoder.append(DecoderBlock(self._conv_block(i, m), self._up_block(m, o)))
        i, o = self.nb_features[0], self.nb_features[0]
        if self.nb_levels > 1 and self.combine == "cat":
            i *= 2
        decoder.append(self._conv_block(i, o))
        self.decoder = nn.Sequential(*decoder)

    def _conv_block(self, i, o=None):
        return ConvGroup(
            self.ndim,
            i,
            o or i,
            kernel_size=self.kernel_size,
            nb_conv=self.nb_conv,
            activation=self.activation,
            norm=self.norm,
            dropout=self.dropout,
            order=self.order,
            residual=self.residual,
        )

    def _down_block(self, i, o=None):
        return ConvBlockDown(
            self.ndim,
            i,
            o or i,
            factor=self.factor,
            kernel_size=1,
            activation=self.activation,
            norm=self.norm,
            dropout=self.dropout,
            order=self.order,
        )

    def _up_block(self, i, o=None):
        return ConvBlockUp(
            self.ndim,
            i,
            o or i,
            factor=self.factor,
            kernel_size=1,
            activation=self.activation,
            norm=self.norm,
            dropout=self.dropout,
            order=self.order,
            combine=self.combine,
        )

    def forward(self, x):
        nb_levels = len(self.encoder)
        if any(s < 2**nb_levels for s in x.shape[2:]):
            raise ValueError(
                f"UNet with {nb_levels} levels requires input spatial dims >= {2**nb_levels}"
            )
        skips = []
        for layer in self.encoder:
            x = layer(x)
            skips.append(x)
        x = skips.pop(-1)
        for n in range(len(self.decoder) - 1):
            x = self.decoder[n].conv(x)
            x = self.decoder[n].up(x, skips.pop(-1))
        x = self.decoder[-1](x)
        return x


class SegNet(nn.Sequential):
    def __init__(
        self,
        ndim,
        in_channels,
        out_channels,
        kernel_size=3,
        activation="Softmax",
        backbone="UNet",
        kwargs_backbone=None,
    ):
        if isinstance(backbone, str):
            backbone_cls = globals()[backbone]
            backbone = backbone_cls(ndim, **(kwargs_backbone or {}))
        act = None
        if (
            activation
            and isinstance(activation, str)
            and activation.lower() == "softmax"
        ):
            act = nn.Softmax(1)
        feat = ConvBlock(
            ndim,
            in_channels,
            backbone.in_channels,
            kernel_size=kernel_size,
            activation=None,
            norm=None,
        )
        pred = ConvBlock(
            ndim,
            backbone.out_channels,
            out_channels,
            kernel_size=1,
            activation=act,
            norm=None,
        )
        super().__init__(feat, backbone, pred)


class SynthSeg(nn.Module):
    def __init__(self, segnet):
        super().__init__()
        self.segnet = segnet

    def forward(self, x):
        return self.segnet(x)


class Model(nn.Module):
    def __init__(
        self,
        ndim: int = 3,
        nb_classes: int = 9,
        seg_nb_levels: int = 6,
        seg_features: Sequence[int] = (16, 24, 32, 48, 64, 96),
        seg_activation: str = "LeakyReLU",
        seg_nb_conv: int = 2,
        seg_norm: Optional[str] = "instance",
        **kwargs,
    ):
        super().__init__()
        backbone = UNet(
            ndim,
            nb_levels=seg_nb_levels,
            nb_features=seg_features,
            nb_conv=seg_nb_conv,
            activation=seg_activation,
            norm=seg_norm,
        )
        segnet = SegNet(ndim, 1, nb_classes, backbone=backbone, activation=None)
        self.network = SynthSeg(segnet)

    def forward(self, x):
        return self.network(x)


# add Real synth patch class and apply after instantiation
class Real(torch.nn.Module):
    def __init__(self, synth):
        super().__init__()
        self.synth = synth

    def forward(self, slab, img, lab):
        mask = (lab >= 1001) & (lab <= 1008)
        new_lab = torch.zeros_like(lab)
        new_lab[mask] = lab[mask] - 1000
        return img, new_lab, img, new_lab


# ==============================================================================
# MAIN SCRIPT LOGIC
# ==============================================================================
def load_nifti(path, device="cpu"):
    """Load NIfTI file using cornucopia LoadTransform - the working method."""
    load_transform = LoadTransform(ndim=3, dtype=torch.float32, device=device)
    tensor = load_transform(path)
    print(f"Loaded image - shape: {tensor.shape}, dtype: {tensor.dtype}")
    print(f"Data range: [{tensor.min():.6f}, {tensor.max():.6f}]")
    return tensor


def save_nifti(pred, output_path, reference_path):
    """Save prediction as NIfTI file - matching working inference.py exactly"""
    # Get affine from reference image if available
    affine = np.eye(4)
    if reference_path and os.path.exists(reference_path):
        try:
            reference_img = nib.load(reference_path)
            affine = reference_img.affine
        except Exception as e:
            print(f"Warning: Could not load affine from {reference_path}: {e}")

    # Convert to numpy and ensure correct data type
    pred_np = pred.cpu().numpy()

    # CRITICAL: Apply the same processing as in the working inference.py
    # Flip over x-axis (axis=0) to maintain same alignment as input image
    if pred_np.shape and pred_np.shape[0] > 0:
        pred_np = np.flip(pred_np, axis=0)

    # Create NIfTI image with the same affine as the reference
    nii_img = nib.Nifti1Image(pred_np.astype(np.uint8), affine)

    # Save as NIfTI
    nib.save(nii_img, output_path)


def main():
    start_total_time = time.time()
    parser = argparse.ArgumentParser(
        description="SynthTopo hippocampus segmentation inference"
    )
    parser.add_argument("input", help="Input NIfTI file (.nii or .nii.gz)")
    parser.add_argument("checkpoint", help="Model checkpoint (.ckpt file)")
    parser.add_argument(
        "-o", "--output", required=True, help="Output segmentation file"
    )
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument(
        "--nb-classes",
        type=int,
        default=9,
        help="Number of output classes (default: %(default)s)",
    )
    parser.add_argument(
        "--seg-nb-levels",
        type=int,
        default=6,
        help="Number of UNet levels (default: %(default)s)",
    )
    parser.add_argument(
        "--seg-features",
        type=int,
        nargs="+",
        default=[16, 24, 32, 48, 64, 96],
        help="UNet feature sizes per level, space-separated (default: %(default)s)",
    )
    args = parser.parse_args()
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    if not os.path.exists(args.checkpoint):
        raise FileNotFoundError(f"Checkpoint file not found: {args.checkpoint}")
    device = torch.device(
        "cuda"
        if (args.device == "auto" and torch.cuda.is_available())
        else (args.device if args.device != "auto" else "cpu")
    )
    print(f"Using device: {device}")
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    print("Loading original model...")
    model = Model(
        nb_classes=args.nb_classes,
        seg_nb_levels=args.seg_nb_levels,
        seg_features=tuple(args.seg_features),
    )
    # patch synth like original inference
    if hasattr(model.network, "synth"):
        model.network.synth = Real(getattr(model.network, "synth", None))
    first_param = next(model.parameters()).clone().detach()
    checkpoint = torch.load(args.checkpoint, map_location="cpu")
    if "state_dict" in checkpoint:
        model.load_state_dict(checkpoint["state_dict"])
        print("Loaded state_dict from checkpoint")
    else:
        model.load_state_dict(checkpoint)
        print("Loaded checkpoint directly")
    if torch.allclose(first_param, next(model.parameters())):
        print("WARNING: parameters unchanged after loading checkpoint")
    model.to(device).eval()
    print("Loading input image...")
    img = LoadTransform(ndim=3, dtype=torch.float32, device=device)(args.input)
    if img.dim() == 4:
        img = img.unsqueeze(0)  # (1,1,D,H,W)
    else:
        raise ValueError(f"Unexpected input shape {img.shape}")
    start_inference_time = time.time()
    with torch.no_grad():
        pred_logits = model(img)
    inference_time = time.time() - start_inference_time
    pred = pred_logits.argmax(dim=1) if pred_logits.shape[1] > 1 else pred_logits
    pred = pred[0]
    print(f"Model inference time: {inference_time:.3f} s")
    save_nifti(pred, args.output, args.input)
    total_time = time.time() - start_total_time
    print("Saved segmentation:", args.output)
    print(f"Total processing time: {total_time:.3f} s")


if __name__ == "__main__":
    main()
