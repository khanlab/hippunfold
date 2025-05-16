import os
import urllib.request
import zipfile


def download_extract(unzip_dir, url):

    outdir = str(unzip_dir)
    os.makedirs(outdir, exist_ok=True)

    zip_path = os.path.join(outdir, "temp.zip")
    urllib.request.urlretrieve("https://" + url, zip_path)

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(outdir)

    os.remove(zip_path)


unzip_dir = snakemake.output.unzip_dir
url = snakemake.params.url

download_extract(unzip_dir, url)
