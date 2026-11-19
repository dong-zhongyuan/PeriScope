import unittest
import numpy as np
from PIL import Image
from pdv3.spatial_images import as_rgb,centered_patch


class SpatialImageTests(unittest.TestCase):
    def test_associated_alpha_is_unpremultiplied(self):
        image=Image.frombytes('La',(1,1),bytes([64,128]))
        self.assertTrue(abs(as_rgb(image).getpixel((0,0))[0]-127)<=1)

    def test_crop_keeps_center_and_exact_pixels_without_resampling(self):
        a=np.zeros((10,12,3),np.uint8);a[:,:,0]=np.arange(12)
        im,meta=centered_patch(a,6,5,4,output=4)
        np.testing.assert_array_equal(np.asarray(im),a[3:7,4:8])
        self.assertEqual(meta['padding'],[0,0,0,0])
        self.assertEqual((meta['center_dx'],meta['center_dy']),(0,0))

    def test_edge_is_reflected_without_recentering(self):
        a=np.arange(5*5*3,dtype=np.uint8).reshape(5,5,3)
        im,meta=centered_patch(a,0,0,4,output=4)
        expected=np.pad(a[:2,:2],((2,0),(2,0),(0,0)),mode='reflect')
        np.testing.assert_array_equal(np.asarray(im),expected)
        self.assertEqual(meta['padding'],[2,2,0,0])
        self.assertEqual((meta['center_dx'],meta['center_dy']),(0,0))

    def test_bad_coordinate_and_nonpositive_size_rejected(self):
        a=np.zeros((4,4,3),np.uint8)
        for x,y,side in [(-1,1,2),(4,1,2),(1,np.nan,2),(1,1,0)]:
            with self.assertRaises(ValueError):centered_patch(a,x,y,side)


if __name__=='__main__':unittest.main()
