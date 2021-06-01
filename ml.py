from fastai.vision.all import *
import torch

path = "images"

fnames = get_image_files(path)

def label_func(x): return x.parent.name

dls = ImageDataLoaders.from_path_func(path, fnames, label_func, valid_pct = 0.5)