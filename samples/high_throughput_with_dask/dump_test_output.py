import pickle
import base64
import sys

base64_data = sys.stdin.read()
pickle_data = base64.b64decode(base64_data)
data = pickle.loads(pickle_data)

print(data)
