# Convert a Blender PNG render to raw RGBA bytes + a dims file, so R (no png package on the host) can place it.
import sys, json
from PIL import Image
im = Image.open(sys.argv[1]).convert('RGBA'); open(sys.argv[2], 'wb').write(im.tobytes())
json.dump({'width': im.width, 'height': im.height}, open(sys.argv[2] + '.json', 'w'))
print('RGBA', im.width, im.height)
