# Note: with SciViews Box 2019, one needs to do first in a terminal:
# sudo pip3 install scikit-image==0.14.2

from skimage import data, util
from skimage.measure import label
img = util.img_as_ubyte(data.coins()) > 110
label_img = label(img, connectivity=img.ndim)
import skimage
props = skimage.measure.regionprops(label_img)
# centroid of first labeled object
props[0].centroid
# centroid of first labeled object
props[0]['centroid']

props[0].area
props[0].equivalent_diameter
