"""
Tests configuration template file.
"""
import os
import json

def test_json():
    """
    Tests correct config JSON syntax.
    """
    with open("tauren_config.json", 'r') as conf:
        config = json.load(conf)
    
    return  
