"""
This module includes unit tests of terrain_data module
"""

import terrain_data


def test_new_terrain_data_object():
    """
    Test of simple creating TerrainData object
    """
    terrain = terrain_data.TerrainData("conf.yaml")
    assert terrain is not None


def test_load_yaml_conf():
    """
    Test loading configuration file from yaml file
    """
    terrain = terrain_data.TerrainData("conf.yaml")
    terrain.load_conf_from_yaml()
    assert terrain.conf.has_key('area') == True and \
        terrain.conf.has_key('terrain') == True and \
        terrain.conf.has_key('rivers') == True


def test_load_terrain_data():
    """
    Simple test loading terrain data
    """
    terrain = terrain_data.TerrainData("conf.yaml")
    terrain.load_conf_from_yaml()
    terrain.load_terrain()
    assert len(terrain.terrain_data) > 0


def test_load_rivers_data():
    """
    Simple test loading rivers data
    """
    terrain = terrain_data.TerrainData("conf.yaml")
    terrain.load_conf_from_yaml()
    terrain.load_rivers()
    assert len(terrain.rivers_data_2d) > 0


def test_load_area_data():
    """
    Simple test loading area data
    """
    terrain = terrain_data.TerrainData("conf.yaml")
    terrain.load_conf_from_yaml()
    terrain.load_area()
    assert len(terrain.area_borders_2d) > 0
