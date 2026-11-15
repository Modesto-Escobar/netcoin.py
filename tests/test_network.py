"""The JSON and the files of the page of a network."""
import json
import os
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from netcoin import netCoin
from netcoin.network import getRawName, netCreate, netWrapper, networkJSON


class TestGetRawName(unittest.TestCase):
    """Tests for getRawName filename encoding function."""

    def test_simple_filename(self):
        result = getRawName("/path/to/image.png")
        self.assertTrue(result.endswith(".png"))
        self.assertNotIn("/", result)
        self.assertNotIn("path", result)

    def test_preserves_extension(self):
        for ext in [".jpg", ".gif", ".svg", ".PNG"]:
            self.assertTrue(getRawName(f"image{ext}").endswith(ext))

    def test_different_filenames_produce_different_hashes(self):
        self.assertNotEqual(getRawName("image1.png"), getRawName("image2.png"))

    def test_encodes_the_name_as_hexadecimal(self):
        self.assertEqual(getRawName("a.b.png"), "612e62.png")


def _net(**kwargs):
    nodes = pd.DataFrame({"name": ["A", "B", "C"], "category": ["X", "Y", "X"]})
    links = pd.DataFrame({"Source": ["A", "B"], "Target": ["B", "C"], "weight": [0.5, 0.8]})
    return netCoin(nodes, links, **kwargs)


class TestNetworkJSON(unittest.TestCase):
    def test_returns_dict_with_expected_keys(self):
        result = networkJSON(_net())
        for key in ("nodes", "nodenames", "links", "linknames", "options"):
            self.assertIn(key, result)

    def test_nodes_and_links_are_lists_of_columns(self):
        result = networkJSON(_net())
        self.assertEqual(result["nodenames"], ["name", "category"])
        self.assertEqual(result["nodes"][0], ["A", "B", "C"])
        self.assertEqual(result["linknames"], ["Source", "Target", "weight"])

    def test_link_ends_are_zero_based_positions(self):
        result = networkJSON(_net())
        self.assertEqual(result["links"][0], [0, 1])
        self.assertEqual(result["links"][1], [1, 2])

    def test_node_names_are_strings(self):
        net = netCoin(pd.DataFrame({"name": [1, 2, 3]}), pd.DataFrame({"Source": [1], "Target": [2]}))
        # 1-based positions here, since the links are numbers
        for name in networkJSON(net)["nodes"][0]:
            self.assertIsInstance(name, str)

    def test_tree_with_two_parents_is_skipped(self):
        net = _net(tree=[["A", "B"], ["C", "B"]])
        with self.assertWarns(UserWarning):
            result = networkJSON(net)
        self.assertNotIn("tree", result)

    def test_valid_tree_is_written_with_positions(self):
        result = networkJSON(_net(tree=[["A", "B"], ["A", "C"]]))
        self.assertEqual(result["tree"], [[0, 0], [1, 2]])

    def test_no_links_means_no_links_key(self):
        result = networkJSON(netCoin(pd.DataFrame({"name": ["A", "B"]})))
        self.assertNotIn("links", result)
        self.assertNotIn("linknames", result)

    def test_layout_planes_are_written(self):
        planes = {"one": np.array([[0, 0], [1, 1], [2, 0]]), "two": np.array([[1, 1], [2, 2], [3, 1]])}
        result = networkJSON(_net(layout=planes))
        self.assertEqual(set(result["layouts"]), {"one", "two"})
        self.assertEqual(result["layouts"]["two"][1], [2, 2])

    def test_statistics_adds_types_and_categories(self):
        nodes = pd.DataFrame({"name": ["A", "B"], "g": pd.Categorical(["u", "v"]), "n": [1.0, 2.0]})
        result = networkJSON(netCoin(nodes, statistics=True))
        self.assertEqual(result["nodetypes"], ["string", "string", "number"])
        self.assertEqual(result["options"]["categories"], {"g": ["u", "v"]})


class TestNetWrapper(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_without_images_returns_valid_data(self):
        result = netWrapper(_net(), self.test_dir)
        self.assertIn("nodes", result)
        self.assertIn("options", result)
        self.assertFalse(os.path.exists(os.path.join(self.test_dir, "images")))

    def _image(self, name="test.png"):
        path = os.path.join(self.test_dir, name)
        with open(path, "w") as f:
            f.write("fake image")
        return path

    def test_images_are_copied_with_a_safe_name(self):
        img = self._image()
        nodes = pd.DataFrame({"name": ["A", "B"], "photo": [img, img]})
        net = netCoin(nodes, pd.DataFrame({"Source": ["A"], "Target": ["B"]}), image="photo")
        out = os.path.join(self.test_dir, "out")
        os.makedirs(out)
        result = netWrapper(net, out)
        self.assertEqual(len(os.listdir(os.path.join(out, "images"))), 1)
        self.assertEqual(result["options"]["imageItems"], ["photo_url"])
        self.assertEqual(result["options"]["imageNames"], ["photo"])
        names = result["nodenames"]
        self.assertEqual(result["nodes"][names.index("photo")], ["test", "test"])
        self.assertTrue(result["nodes"][names.index("photo_url")][0].startswith("images/"))

    def test_several_images_in_one_cell(self):
        a, b = self._image("a.png"), self._image("b.png")
        nodes = pd.DataFrame({"name": ["A"], "photo": [a + "|" + b]})
        net = netCoin(nodes, image="photo")
        out = os.path.join(self.test_dir, "out")
        os.makedirs(out)
        result = netWrapper(net, out)
        self.assertEqual(len(os.listdir(os.path.join(out, "images"))), 2)
        url = result["nodes"][result["nodenames"].index("photo_url")][0]
        self.assertEqual(url.count("|"), 1)

    def test_a_missing_image_warns(self):
        nodes = pd.DataFrame({"name": ["A"], "photo": ["/no/such/file.png"]})
        out = os.path.join(self.test_dir, "out")
        os.makedirs(out)
        with self.assertWarns(UserWarning):
            netWrapper(netCoin(nodes, image="photo"), out)

    def test_the_network_is_not_modified(self):
        img = self._image()
        nodes = pd.DataFrame({"name": ["A"], "photo": [img]})
        net = netCoin(nodes, image="photo")
        out = os.path.join(self.test_dir, "out")
        os.makedirs(out)
        netWrapper(net, out)
        self.assertEqual(net.nodes["photo"].tolist(), [img])
        self.assertEqual(net.options["imageItems"], ["photo"])
        self.assertNotIn("imageNames", net.options)

    def test_background_image_is_copied(self):
        img = self._image("bg.png")
        net = netCoin(pd.DataFrame({"name": ["A"]}), background=img)
        out = os.path.join(self.test_dir, "out")
        os.makedirs(out)
        result = netWrapper(net, out)
        self.assertTrue(result["options"]["background"].startswith('url("images/'))
        self.assertEqual(len(os.listdir(os.path.join(out, "images"))), 1)

    def test_statistics_icons_are_copied(self):
        out = os.path.join(self.test_dir, "out")
        os.makedirs(out)
        netWrapper(netCoin(pd.DataFrame({"name": ["A"]}), statistics=True), out)
        self.assertTrue(os.path.isfile(os.path.join(out, "stats-icons", "anova.png")))


class TestNetworkIntegration(unittest.TestCase):
    def test_the_page_is_written(self):
        folder = tempfile.mkdtemp()
        try:
            out = os.path.join(folder, "net")
            netCreate(_net(main="Test Network"), out)
            html = open(os.path.join(out, "index.html"), encoding="utf-8").read()
            data = json.loads(html.split('id="data">')[1].split("</script>")[0])
            self.assertEqual(data["options"]["main"], "Test Network")
            self.assertEqual(len(data["nodes"]), 2)
            self.assertEqual(len(data["nodes"][0]), 3)
        finally:
            shutil.rmtree(folder, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
