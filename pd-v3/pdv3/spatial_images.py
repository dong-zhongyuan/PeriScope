"""Verified TIFF streaming and coordinate-preserving spatial image patches."""
from contextlib import contextmanager
import hashlib
from pathlib import Path
import tempfile

import numpy as np
from PIL import Image


def as_rgb(image):
    # Pillow requires an explicit unpremultiply step for associated-alpha TIFFs.
    # Ignoring the alpha convention changes intensities at transparent edges.
    if image.mode=='La':image=image.convert('LA')
    elif image.mode=='RGBa':image=image.convert('RGBA')
    return image.convert('RGB')


@contextmanager
def tiff_member(archive, info, scratch):
    """Use one temporary TIFF at a time; zipfile checks its CRC at EOF."""
    scratch=Path(scratch);scratch.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='tiff_',dir=scratch) as tmp:
        path=Path(tmp)/'image.tif';sha=hashlib.sha256();n=0
        with archive.open(info) as inp,path.open('wb') as out:
            while chunk:=inp.read(8*1024*1024):out.write(chunk);sha.update(chunk);n+=len(chunk)
        if n!=info.file_size:raise ValueError('Incomplete TIFF')
        yield path,dict(member=info.filename,bytes=n,sha256=sha.hexdigest(),crc32=f'{info.CRC:08x}',zip_crc_verified_to_eof=True)


def centered_patch(rgb, x, y, side, output=224):
    """Keep the requested center; reflect only pixels beyond image edges.

    Integer crop placement differs by at most half a pixel from the supplied
    center. Spatial selection never depends on RNA or disease labels.
    """
    if rgb.ndim!=3 or rgb.shape[2]!=3 or rgb.dtype!=np.uint8:raise ValueError('Expected uint8 RGB image')
    height,width=rgb.shape[:2]
    if not np.isfinite([x,y,side]).all() or not (0<=x<width and 0<=y<height) or side<=0:raise ValueError('Invalid patch geometry')
    length=max(1,int(round(side)));left=int(round(x-length/2));top=int(round(y-length/2))
    right=left+length;bottom=top+length
    pad=(max(0,-left),max(0,-top),max(0,right-width),max(0,bottom-height))
    crop=rgb[max(0,top):min(height,bottom),max(0,left):min(width,right)]
    if any(pad):crop=np.pad(crop,((pad[1],pad[3]),(pad[0],pad[2]),(0,0)),mode='reflect')
    if crop.shape[:2]!=(length,length):raise ValueError('Patch geometry was not preserved')
    image=Image.fromarray(crop).resize((output,output),Image.Resampling.BICUBIC)
    return image,dict(left=left,top=top,side=length,padding=list(pad),center_dx=left+length/2-x,center_dy=top+length/2-y)
