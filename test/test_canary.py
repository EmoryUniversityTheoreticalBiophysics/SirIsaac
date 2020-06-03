# 1.16.2019
#
# Bryan Daniels
#
# Simplest test that the package loads, etc
#

import unittest

import SirIsaac

class TestCanary(unittest.TestCase):
    def test_addition(self):
        """
        Test 1 + 1 = 2
        """
        self.assertEqual(2, 1+1)
        
