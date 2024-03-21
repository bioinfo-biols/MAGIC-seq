import math
import numpy as np
def load_image(filename):
    print(f'Image loading from {filename}')
    from PIL import Image
    from PIL import ImageFile
    ImageFile.LOAD_TRUNCATED_IMAGES = True
    Image.MAX_IMAGE_PIXELS = None
    img = Image.open(filename)
    img = np.array(img)
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
        # remove alpha channel
    return img

def hex2rgb(hex_str):
    h = hex_str.lstrip('#')
    rgb = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    r,g,b=rgb
    rgb=(r/255,g/255,b/255)
    return rgb

#return grey image
def rgb2gray(rgb):
    import numpy as np
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140]).astype('uint8') 
