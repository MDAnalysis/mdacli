# Adapted from mdtraj/devtools/travis-ci/update_versions_json.py,
# by Robert T McGibbon, under the LGPLv2.1 license
# Similarly by the terms of LGPL this is also in the public domain
# Lily Wang, 2020
#
# This is a hacky DIY of ReadTheDoc's versioned docs
# First we use msmb_theme for the little versions up
# This relies on a `versions.json` at the root or `base_url`
# e.g. mdacli.mdanalysis.org/versions.json
# versions.json consists of a List[Dict[str, str]] with the following keys:
# version: version label
# display: how to display in version popup
# url: URL to documentation for that version (a subfolder)
# latest: whether it's the latest version

# The docs are pushed to the gh-pages branch
# In the CI we make a version subfolder and
# move the built documentation into the folder.
# In this Python script we then update versions.json
# with the version URL. We also determine if this is
# the latest or development version.

# If so, we copy the docs to a new subfolder:
# stable/ or dev/ . That's so we have a stable endpoint
# for the most up-to-date documentation at any point.
# For stable, we also write root-level redirects.
# i.e. mdacli.mdanalysis.org/index.html will redirect
#      to mdacli.mdanalysis.org/stable/index.html

# Finally, we create a sitemap index of
# all the individual version sitemaps.


import json
import os
import shutil
import xml.etree.ElementTree as ET
import errno
import glob
import textwrap
import shutil
import argparse
from typing import Callable, Optional, Any

try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen

parser = argparse.ArgumentParser("Update JSON stubs and sitemaps for docs")
parser.add_argument("--url", required=True, type=str)
parser.add_argument("--version", required=True, type=str)
parser.add_argument("--devlabel", default="dev", required=False,
                    help="marker to show development version")


class VersionJazz:
    """Class to manage version stuff"""

    def __init__(self, base_url: str, version: str, devlabel: str = "dev"):
        self.base_url = base_url
        self.version = version
        self.devlabel = devlabel
        self.is_latest = devlabel not in version

        self.versions = self.get_web_file(filename="versions.json",
                                          callback=json.loads,
                                          default=[])
        self.latest_version = self.get_latest_version()
        self.development_version = self.get_development_version()

    def update_versions(self):
        """Update versions.json with new version
        and add HTML redirects to version directories

        The steps followed here are:

        1. Add current version to ``self.versions``
        2. Sort versions
        3. Add redirect stubs to latest and development versions
        4. Write file out to versions.json
        """
        # add new version
        self.add_current_version()
        self.latest_version = self.get_latest_version()
        self.development_version = self.get_development_version()
        self.versions.sort(key=lambda x: x["version"])

        self.add_redirect_stubs_to_stable_and_dev()
        self.dump_versions()

    def get_web_file(
            self,
            filename: str,
            callback: Callable = json.loads,
            default: Any = None,
            ) -> Any:
        """Get file at URL, applying a callback to it to return data

        If an error is raised, the ``default`` is returned.
        The URL is created from ``self.base_url`` and ``filename``.
        """
        url = os.path.join(self.base_url, filename)
        try:
            page = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            data = urlopen(page).read().decode()
        except Exception as e:
            print(e)
            try:
                with open(filename, 'r') as f:
                    return callback(f)
            except IOError as e:
                print(e)
                return default
        else:
            return callback(data)

    def write_redirect(
            self,
            source_version: str,
            source_file: str = "",
            target: Optional[str] = None,
            ):
        """Write HTML redirect from source_version/source_file to target"""

        if target is None:
            target = source_file
        url = os.path.join(self.base_url, source_version, source_file)
        REDIRECT = textwrap.dedent(f"""
        <!DOCTYPE html>
        <meta charset="utf-8">
        <title>Redirecting to {url}</title>
        <meta http-equiv="refresh" content="0; URL={url}">
        <link rel="canonical" href="{url}">
        """)
        with open(target, 'w') as f:
            f.write(REDIRECT)
        print(f"Wrote redirect from {url} to {target}")

    def get_version_by_criteria(self, callable: Callable) -> Optional[str]:
        """Find version satisfying criteria, or return last version"""
        for ver in self.versions[::-1]:
            if callable(ver):
                return ver["version"]
        try:
            return self.versions[-1]["version"]
        except IndexError:
            return None

    def get_latest_version(self) -> Optional[str]:
        return self.get_version_by_criteria(lambda x: x["latest"])

    def get_development_version(self) -> Optional[str]:
        def is_dev(x):
            return (self.devlabel in x["version"]
                    and not x["version"] == self.devlabel)
        return self.get_version_by_criteria(is_dev)

    def add_current_version(self):
        """Update version dictionary **in place**"""
        existing = [item["version"] for item in self.versions]

        if self.version not in existing:
            if self.is_latest:
                # mark old versions *not* latest
                for ver in self.versions:
                    ver["latest"] = False
            self.versions.append({
                "version": self.version,
                "display": self.version,
                "url": os.path.join(self.base_url, self.version),
                "latest": self.is_latest
                })

    def redirect_sitemap(self, target: str = "stable"):
        """Replace paths in copied sitemap.xml with new directory path

        Sitemaps can only contain URLs 'within' that directory structure.
        For more, see https://www.sitemaps.org/faq.html#faq_sitemap_location
        """

        sitemap = os.path.join(target, "sitemap.xml")
        old_url = os.path.join(self.base_url, self.version)
        new_url = os.path.join(self.base_url, target)

        try:
            with open(sitemap, "r") as f:
                contents = f.read()
        except OSError:
            raise ValueError(f"{sitemap} not found")
        redirected = contents.replace(old_url, new_url)
        with open(sitemap, "w") as f:
            f.write(redirected)
        print(f"Redirected URLs in {sitemap} from {old_url} to {new_url}")

    def add_or_update_version(self, version: str):
        """Add or update the version path to self.versions **in-place**"""
        for ver in self.versions:
            if ver["version"] == version:
                ver["url"] = os.path.join(self.base_url, version)
                break
        else:
            self.versions.append({
                "version": version,
                "display": version,
                "url": os.path.join(self.base_url, version),
                "latest": False
                })

    def copy_version(self, target: str = "stable"):
        """Copy a version directory to target and update sitemap URLs"""

        shutil.copytree(self.version, target)
        print(f"Copied {self.version} to {target}")
        self.redirect_sitemap(target)
        self.add_or_update_version(target)

    def write_root_redirects_to_target(self, target: str = "stable"):
        """Redirect base.url/file to base.url/target/file"""

        html_files = glob.glob(f"{target}/**/*.html", recursive=True)
        for file in html_files:
            # below should be true because we only globbed stable/* paths
            assert file.startswith(target)
            root_url = file[len(target):]
            while root_url[0] == "/":
                root_url = root_url[1:]

            # check directory exists
            dirname = os.path.dirname(root_url)
            if dirname and not os.path.exists(dirname):
                try:
                    os.makedirs(dirname)
                except OSError as exc:
                    if exc.errno != errno.EEXIST:
                        raise

            self.write_redirect(
                source_file=file, source_version="", target=root_url)

    def add_redirect_stubs_to_stable_and_dev(self):
        """Add the stable, latest, development docs by copying"""

        if self.is_latest:
            self.copy_version("stable")
            self.write_root_redirects_to_target("stable")

        if self.latest_version:
            self.write_redirect(source_version=self.latest_version,
                                source_file="index.html",
                                target="latest/index.html")

        if self.development_version == self.version:
            self.copy_version("dev")

        self.write_redirect(source_version="stable",
                            source_file="index.html")

    def dump_versions(self):
        """Write versions out again"""
        with open("versions.json", "w") as f:
            json.dump(self.versions, f, indent=2)

    def write_overall_sitemap(self, filename: str = "sitemap_index.xml"):
        """Write sitemap index for individual version sitemaps"""
        ET.register_namespace('xhtml', "http://www.w3.org/1999/xhtml")

        # so we could make 1 big sitemap as commented
        # below, but they must be max 50 MB / 50k URL.
        # Yes, this is 100+ releases, but who knows when
        # that'll happen and who'll look at this then?
        # bigroot = ET.Element("urlset")
        # bigroot.set("xmlns", "http://www.sitemaps.org/schemas/sitemap/0.9")
        # for ver in versions:
        #     tree = get_web_file(ver['version']+'/sitemap.xml', ET.fromstring,
        #                         ET.fromstring(''))
        #     root = tree.getroot()
        #     bigroot.extend(root.getchildren())
        # ET.ElementTree(bigroot).write('sitemap.xml',
        #                               xml_declaration=True,
        #                               encoding='utf-8',
        #                               method="xml")

        # so instead we make a sitemap of sitemaps.
        bigroot = ET.Element("sitemapindex")
        bigroot.set("xmlns", "http://www.sitemaps.org/schemas/sitemap/0.9")
        for ver in self.versions:
            path = os.path.join(self.base_url, f"{ver['version']}/sitemap.xml")
            sitemap = ET.SubElement(bigroot, "sitemap")
            ET.SubElement(sitemap, "loc").text = path

        ET.ElementTree(bigroot).write(filename,
                                      xml_declaration=True,
                                      encoding="utf-8",
                                      method="xml")
        print(f"Wrote to {filename}")


if __name__ == "__main__":
    args = parser.parse_args()

    # check values
    if "http" not in args.url:
        raise ValueError("URL should have the transfer protocol (HTTP/S). "
                         f"Given: --url {args.url}")
    try:
        int(args.version[0])
    except ValueError:
        raise ValueError("$VERSION should start with a number. "
                         f"Given: --version {args.version}") from None

    jazzer = VersionJazz(base_url=args.url,
                         version=args.version,
                         devlabel=args.devlabel)
    jazzer.update_versions()
    jazzer.write_overall_sitemap()
