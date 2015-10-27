"""
This module includes unit tests of terrain_data module
"""

import terrain_data

def test_new_terrain_data():
	"""
	Test of simple creating TerrainData object
	"""
	terrain = terrain_data.TerrainData("conf.yaml")
	assert terrain != None
