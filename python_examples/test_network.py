import unittest
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from netcoin.network import networkJSON, getRawName, netWrapper, netCreate


class TestGetRawName(unittest.TestCase):
    """Tests for getRawName filename encoding function."""

    def test_simple_filename(self):
        """Should encode filename to hex and preserve extension."""
        result = getRawName("/path/to/image.png")
        self.assertTrue(result.endswith(".png"))
        self.assertNotIn("/", result)
        self.assertNotIn("path", result)

    def test_preserves_extension(self):
        """Should preserve original file extension."""
        for ext in [".jpg", ".gif", ".svg", ".PNG"]:
            filepath = f"image{ext}"
            result = getRawName(filepath)
            self.assertTrue(result.endswith(ext))

    def test_different_filenames_produce_different_hashes(self):
        """Different filenames should produce different hex values."""
        name1 = getRawName("image1.png")
        name2 = getRawName("image2.png")
        self.assertNotEqual(name1, name2)


class TestNetworkJSON(unittest.TestCase):
    """Tests for networkJSON data transformation function."""

    def setUp(self):
        """Create sample network data for testing."""
        self.net = {
            'links': {'Source': [0, 1], 'Target': [1, 2], 'weight': [0.5, 0.8]},
            'nodes': {'name': ['A', 'B', 'C'], 'category': ['X', 'Y', 'X']},
            'options': {'nodeName': 'name'},
            'nodeAttrNames': ['name', 'category'],
            'linkAttrNames': ['Source', 'Target', 'weight']
        }

    def test_returns_dict_with_expected_keys(self):
        """Should return data dict with required keys."""
        result = networkJSON(self.net)
        self.assertIn('nodes', result)
        self.assertIn('nodenames', result)
        self.assertIn('links', result)
        self.assertIn('linknames', result)
        self.assertIn('options', result)

    def test_node_names_converted_to_strings(self):
        """Should convert node names to strings."""
        self.net['nodes']['name'] = [1, 2, 3]  # integers
        result = networkJSON(self.net)
        # Check that node names are strings
        for name in result['nodes'][0]:  # First node attribute (names)
            self.assertIsInstance(name, str)

    def test_invalid_tree_is_skipped(self):
        """Should skip tree if it has duplicate parents."""
        self.net['tree'] = [['A', 'B'], ['A', 'B']]  # Duplicate edges
        result = networkJSON(self.net)
        # Should not have tree key if invalid
        self.assertNotIn('tree', result)

    def test_links_none_excluded(self):
        """Should exclude links from output if None."""
        self.net['links'] = None
        result = networkJSON(self.net)
        self.assertNotIn('links', result)


class TestNetWrapper(unittest.TestCase):
    """Tests for netWrapper image and path handling."""

    def setUp(self):
        """Create temporary directory for test files."""
        self.test_dir = tempfile.mkdtemp()
        self.net_base = {
            'links': {'Source': [0], 'Target': [1]},
            'nodes': {'name': ['A', 'B']},
            'options': {'nodeName': 'name'},
            'nodeAttrNames': ['name'],
            'linkAttrNames': ['Source', 'Target']
        }

    def tearDown(self):
        """Clean up test directory."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_without_images_returns_valid_data(self):
        """Should return valid data when no images specified."""
        result = netWrapper(self.net_base, self.test_dir)
        self.assertIn('nodes', result)
        self.assertIn('options', result)

    def test_creates_images_directory_when_needed(self):
        """Should create images directory when images are present."""
        # Create a temporary image file
        img_path = os.path.join(self.test_dir, "test.png")
        with open(img_path, 'w') as f:
            f.write("fake image")

        net = self.net_base.copy()
        net['nodes']['photo'] = [img_path]  # Column with file paths
        net['options']['imageItems'] = ['photo']  # Column name, NOT the file path
        net['options']['imageNames'] = ['photo']  # Original column name

        result = netWrapper(net, self.test_dir)
        img_dir = os.path.join(self.test_dir, 'images')
        self.assertTrue(os.path.exists(img_dir))

        # Verify image was copied
        self.assertGreater(len(os.listdir(img_dir)), 0)


class TestNetworkIntegration(unittest.TestCase):
    """Integration tests for complete network creation workflow."""

    def setUp(self):
        """Create test data."""
        self.test_dir = tempfile.mkdtemp()
        self.net = {
            'links': {'Source': [0, 1], 'Target': [1, 2]},
            'nodes': {'name': ['A', 'B', 'C'], 'size': [10, 20, 15]},
            'options': {'nodeName': 'name', 'main': 'Test Network'},
            'nodeAttrNames': ['name', 'size'],
            'linkAttrNames': ['Source', 'Target']
        }

    def tearDown(self):
        """Clean up."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_full_pipeline_produces_valid_output(self):
        """Full wrap -> JSON pipeline should produce valid structure."""
        result = netWrapper(self.net, self.test_dir)

        # Verify structure
        self.assertEqual(len(result['nodes']), 2)  # name + size
        self.assertEqual(len(result['nodenames']), 2)
        self.assertEqual(result['options']['main'], 'Test Network')

        # Verify node data is properly organized
        names = result['nodes'][0]
        sizes = result['nodes'][1]
        self.assertEqual(len(names), 3)
        self.assertEqual(len(sizes), 3)


if __name__ == '__main__':
    unittest.main()
