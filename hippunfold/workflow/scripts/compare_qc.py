from PIL import Image
import os

ref_path = (
    "ref_data/qc/sub-01_space-cropT1w_desc-subfields_atlas-multihist7_volumes.png"
)
new_path = "/mnt/test_out/sub-01/qc/sub-01_space-cropT1w_desc-subfields_atlas-multihist7_volumes.png"
out_path = "sub-01_qc_comparison.png"

if not os.path.exists(ref_path):
    raise FileNotFoundError(f"Reference image not found: {ref_path}")
if not os.path.exists(new_path):
    raise FileNotFoundError(f"New image not found: {new_path}")

ref = Image.open(ref_path)
new = Image.open(new_path)

combined = Image.new("RGB", (ref.width + new.width, max(ref.height, new.height)))
combined.paste(ref, (0, 0))
combined.paste(new, (ref.width, 0))

combined.save(out_path)
print(f"Saved combined QC image to {out_path}")
