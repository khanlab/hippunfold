# populate the HIPPUNFOLD_CACHE_DIR folder as needed
from lib import utils as utils
from snakebids_templates import generate_templates

download_dir = utils.get_download_dir()


# zenodo.org URLs are handled by the cached-http plugin (rate limits, retries,
# checksums); keep its cache inside HIPPUNFOLD_CACHE_DIR instead of ~/.cache
storage cached_http:
    provider="cached-http",
    cache=str(Path(download_dir) / "http_cache"),


storage http:
    provider="http",


def template_storage(url):
    """cached-http only handles zenodo.org, so fetch other URLs (e.g. OSF) with http"""
    if "zenodo.org" in url:
        return storage.cached_http(url)
    return storage.http(url)


# 3D templates, 3D atlases, 2D templates and 2D atlases (see templates in the config)
templates = generate_templates(config["templates"], workflow.basedir)

COMPONENTS_3D = [c for c in templates if c.startswith(("template3d", "atlas3d"))]
COMPONENTS_2D = [c for c in templates if c.startswith(("template2d", "atlas2d"))]

# 3D files are imported here; 2D files into the atlas dir, alongside atlases created
# with group_create_atlas
template_dir = Path(download_dir) / "templates"


# run on the host (not submitted to a cluster), since compute nodes may lack internet access
localrules:
    import_template_file,
    import_atlas_file,


rule import_template_file:
    """copy a 3D template or atlas file into the cache dir, with its canonical name"""
    input:
        templates.import_input(template_storage, COMPONENTS_3D),
    output:
        template_dir / "{relpath}",
    wildcard_constraints:
        relpath=templates.import_constraint(COMPONENTS_3D),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(10),
    run:
        templates.import_file(input[0], output[0], wildcards.relpath)


rule import_atlas_file:
    """extract a 2D template or atlas file into the atlas dir"""
    input:
        templates.import_input(template_storage, COMPONENTS_2D),
    output:
        get_atlas_dir() / "{relpath}",
    wildcard_constraints:
        relpath=templates.import_constraint(COMPONENTS_2D),
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(10),
    run:
        templates.import_file(input[0], output[0], wildcards.relpath)


def template_file(component, **entities):
    """path of an imported 3D template file with the given (fixed) entities"""
    return str(template_dir / templates[component].relpath(**entities))


def has_template_file(component, **entities):
    return len(templates[component].filter(**entities)) > 0


def _template_hemis(component, wildcards, entities):
    """hemispheres of the template file matching the entities and rule wildcards

    For files without a hemi entity (e.g. whole-head images), these are the
    hemispheres of the template's 3D atlas."""
    comp = templates[component]
    if "hemi" not in comp.wildcards:
        return templates["atlas3d_dseg"].filter(tpl=entities["tpl"]).entities["hemi"]
    query = {
        w: entities[w] if w in entities else wildcards[w]
        for w in comp.wildcards
        if w != "hemi" and (w in entities or w in wildcards.keys())
    }
    return comp.filter(**query).entities.get("hemi", [])


def template_hemi(component, **entities):
    """the rule's hemi if the template has it, otherwise the other hemi (to flip)

    Templates may have only a single hemisphere, or be missing files for one."""

    def get(wildcards):
        hemis = _template_hemis(component, wildcards, entities)
        if wildcards.hemi in hemis:
            return wildcards.hemi
        if hemis:
            return hemis[0]
        raise ValueError(f"No {component} file in the template for {entities}")

    return get


def imported_template_file(component, **entities):
    """input function for an imported 3D template file, flipped from the other
    hemisphere if needed"""
    if "hemi" in templates[component].wildcards:
        entities.setdefault("hemi", template_hemi(component, **entities))
    return templates[component].imported(template_dir, **entities)


def copy_or_flip(component, **entities):
    """params function: copy the template file, or flip it if it is from the other hemi"""

    def get(wildcards, input, output):
        if wildcards.hemi in _template_hemis(component, wildcards, entities):
            return f"cp {input} {output}"
        return f"c3d {input} -flip x -o {output}"

    return get


## unpack template
# this is used in both shape_inject.smk and templateseg.smk, so instead of having redundant rules I will simply apply them hemisphere
# Template-based segmentation supports templates that have only a single hemisphere
# by flipping it


def crop_ref_suffix(template):
    """contrast of the template's corobl crop reference: T2w if available"""
    if has_template_file("template3d_crop", tpl=template, suffix="T2w"):
        return "T2w"
    return "T1w"


rule import_template_dseg:
    input:
        imported_template_file(
            "atlas3d_dseg", tpl=config["inject_template"], desc="hipptissue"
        ),
    output:
        template_seg=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                desc="hipptissue",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
            )
        ),
    group:
        "shapeinject"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        cmd=copy_or_flip(
            "atlas3d_dseg", tpl=config["inject_template"], desc="hipptissue"
        ),
    shell:
        "{params.cmd}"


rule import_template_dseg_dentate:
    input:
        imported_template_file(
            "atlas3d_dseg", tpl=config["inject_template"], desc="dentatetissue"
        ),
    output:
        template_seg=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                desc="dentatetissue",
                hemi="{hemi}",
                suffix="dseg.nii.gz",
            )
        ),
    group:
        "shapeinject"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        cmd=copy_or_flip(
            "atlas3d_dseg", tpl=config["inject_template"], desc="dentatetissue"
        ),
    shell:
        "{params.cmd}"


rule import_template_coords:
    input:
        imported_template_file("atlas3d_coords", tpl=config["inject_template"]),
    output:
        template_coords=temp(
            bids(
                root=root,
                datatype="coords",
                **inputs.subj_wildcards,
                dir="{dir}",
                label="{label}",
                suffix="coords.nii.gz",
                desc="init",
                space="template",
                hemi="{hemi}",
            )
        ),
    group:
        "shapeinject"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        cmd=copy_or_flip("atlas3d_coords", tpl=config["inject_template"]),
    shell:
        "{params.cmd}"


rule import_template_anat:
    input:
        imported_template_file(
            "template3d",
            tpl=config["inject_template"],
            suffix=get_modality_suffix(config["modality"]),
        ),
    output:
        template_anat=temp(
            bids(
                root=root,
                datatype="anat",
                space="template",
                **inputs.subj_wildcards,
                hemi="{hemi}",
                suffix="{modality}.nii.gz".format(
                    modality=get_modality_suffix(config["modality"])
                ),
            ),
        ),
    group:
        "shapeinject"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        cmd=copy_or_flip(
            "template3d",
            tpl=config["inject_template"],
            suffix=get_modality_suffix(config["modality"]),
        ),
    shell:
        "{params.cmd}"


rule import_template_anat_crop:  # used only in templateseg workflow
    input:
        imported_template_file(
            "template3d_crop",
            tpl=config["inject_template"],
            suffix=crop_ref_suffix(config["inject_template"]),
        ),
    output:
        template_anat=temp(
            bids(
                root=root,
                datatype="anat",
                desc="template",
                space="corobl",
                **inputs.subj_wildcards,
                hemi="{hemi}",
                suffix="{modality}.nii.gz".format(
                    modality=get_modality_suffix(config["modality"])
                ),
            ),
        ),
    group:
        "shapeinject"
    resources:
        mem_mb=scale_by_attempt(1000),
        runtime=scale_by_attempt(1),
    params:
        cmd=copy_or_flip(
            "template3d_crop",
            tpl=config["inject_template"],
            suffix=crop_ref_suffix(config["inject_template"]),
        ),
    shell:
        "{params.cmd}"
