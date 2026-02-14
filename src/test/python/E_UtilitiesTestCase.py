from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path("../../main/python").resolve()))

from E_Utilities import find_names_or_none

from bs4 import BeautifulSoup


class E_UtilitiesTestCase(unittest.TestCase):
    """Pure unit tests for E_Utilities functions."""

    # -- find_names_or_none tests --

    def test_find_names_or_none_nested_tags(self):
        """Traverses nested tags to find text."""
        xml = "<root><AuthorList><Author><LastName>Smith</LastName></Author></AuthorList></root>"
        soup = BeautifulSoup(xml, "xml")
        result = find_names_or_none(soup, ["AuthorList", "Author", "LastName"])
        self.assertEqual(result, "Smith")

    def test_find_names_or_none_missing_tag(self):
        """Returns None when intermediate tag is missing."""
        xml = "<root><AuthorList><Author><LastName>Smith</LastName></Author></AuthorList></root>"
        soup = BeautifulSoup(xml, "xml")
        result = find_names_or_none(soup, ["AuthorList", "Missing"])
        self.assertIsNone(result)

    def test_find_names_or_none_with_attribute(self):
        """Extracts tag attribute value."""
        xml = '<root><Type value="protein-coding"/></root>'
        soup = BeautifulSoup(xml, "xml")
        result = find_names_or_none(soup, ["Type"], attribute="value")
        self.assertEqual(result, "protein-coding")

    def test_find_names_or_none_single_tag(self):
        """Finds text in a single tag name."""
        xml = "<root><Title>Some Title</Title></root>"
        soup = BeautifulSoup(xml, "xml")
        result = find_names_or_none(soup, ["Title"])
        self.assertEqual(result, "Some Title")

    def test_find_names_or_none_missing_first_tag(self):
        """Returns None when first tag is not found."""
        xml = "<root><Other>text</Other></root>"
        soup = BeautifulSoup(xml, "xml")
        result = find_names_or_none(soup, ["Missing"])
        self.assertIsNone(result)

    def test_find_names_or_none_missing_attribute(self):
        """Returns None when attribute does not exist on found tag."""
        xml = "<root><Type>protein-coding</Type></root>"
        soup = BeautifulSoup(xml, "xml")
        result = find_names_or_none(soup, ["Type"], attribute="nonexistent")
        self.assertIsNone(result)
