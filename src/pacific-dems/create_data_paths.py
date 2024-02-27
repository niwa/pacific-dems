# -*- coding: utf-8 -*-
"""
Created on Friday November 10 11:12:20 2023

Contains routines for creating data paths associated with a country.

@author: pearsonra
"""

import pathlib

def get_paths(country_name: str, resolution: int) -> dict:
    data_path = pathlib.Path().cwd() / ".." / "data"
    fabdem_path = data_path / "FABDEM"
    country_path = data_path / country_name
    country_path.mkdir(parents=True, exist_ok=True)
    output_path = country_path / f"{resolution}m_dems"
    output_path.mkdir(parents=True, exist_ok=True)
    lidar_path = country_path / "lidar"
    land_url = r"https://github.com/digitalearthpacific/depal/raw/main/padm.gpkg"
    
    paths = {"fabdem": fabdem_path, "country": country_path, "output": output_path, "lidar": lidar_path, "land": land_url}
    
    return paths