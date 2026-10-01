"""Read display meshes: binary glTF (.glb), the format NASA publishes in.

Only what a display needs is read -- positions, normals, the triangles, and
each primitive's material name and base colour -- with every node transform
applied, so the result is in the file's own world frame.  Nothing here touches
the physics: a mesh read by this module is drawn, never solved.
"""
from __future__ import annotations

import json
import struct
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

_COMPONENT = {5120: np.int8, 5121: np.uint8, 5122: np.int16, 5123: np.uint16,
              5125: np.uint32, 5126: np.float32}
_WIDTH = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4, "MAT4": 16}


@dataclass
class Primitive:
    """One triangle list from the file, in world coordinates."""

    material: str
    colour: Tuple[float, float, float]
    vertices: np.ndarray             # (n, 3)
    faces: np.ndarray                # (m, 3) indices into vertices
    normals: Optional[np.ndarray]    # (n, 3) or None


def _chunks(path: str):
    with open(path, "rb") as fh:
        data = fh.read()
    magic, version, length = struct.unpack_from("<4sII", data, 0)
    if magic != b"glTF":
        raise ValueError(f"{path} is not a binary glTF (.glb) file")
    if version != 2:
        raise ValueError(f"{path} is glTF version {version}; only 2 is read")
    doc, blob, offset = None, b"", 12
    while offset < length:
        size, kind = struct.unpack_from("<II", data, offset)
        chunk = data[offset + 8: offset + 8 + size]
        if kind == 0x4E4F534A:          # "JSON"
            doc = json.loads(chunk)
        elif kind == 0x004E4942:        # "BIN\0"
            blob = chunk
        offset += 8 + size
    if doc is None:
        raise ValueError(f"{path} has no JSON chunk")
    return doc, blob


def _accessor(doc: dict, blob: bytes, index: int) -> np.ndarray:
    acc = doc["accessors"][index]
    view = doc["bufferViews"][acc["bufferView"]]
    dtype = np.dtype(_COMPONENT[acc["componentType"]])
    width = _WIDTH[acc["type"]]
    start = view.get("byteOffset", 0) + acc.get("byteOffset", 0)
    count = acc["count"]
    item = dtype.itemsize * width
    stride = view.get("byteStride", 0) or item
    rows = np.frombuffer(blob, np.uint8, count=stride * (count - 1) + item, offset=start)
    rows = np.lib.stride_tricks.as_strided(rows, (count, item), (stride, 1))
    out = np.frombuffer(np.ascontiguousarray(rows).tobytes(), dtype).reshape(count, width)
    return out.astype(np.float64 if dtype.kind == "f" else np.int64)


def _local(node: dict) -> np.ndarray:
    if "matrix" in node:
        return np.array(node["matrix"], float).reshape(4, 4).T   # column-major
    t = np.array(node.get("translation", [0.0, 0.0, 0.0]))
    x, y, z, w = node.get("rotation", [0.0, 0.0, 0.0, 1.0])
    s = np.array(node.get("scale", [1.0, 1.0, 1.0]))
    r = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    m = np.eye(4)
    m[:3, :3] = r * s
    m[:3, 3] = t
    return m


def read_glb(path: str) -> List[Primitive]:
    """Every triangle primitive in the default scene, transforms applied."""
    doc, blob = _chunks(path)
    materials = doc.get("materials", [])
    out: List[Primitive] = []

    def walk(index: int, parent: np.ndarray) -> None:
        node = doc["nodes"][index]
        world = parent @ _local(node)
        if "mesh" in node:
            for prim in doc["meshes"][node["mesh"]]["primitives"]:
                if prim.get("mode", 4) != 4:          # triangles only
                    continue
                pos = _accessor(doc, blob, prim["attributes"]["POSITION"])
                verts = (np.c_[pos, np.ones(len(pos))] @ world.T)[:, :3]
                normals = None
                if "NORMAL" in prim["attributes"]:
                    n = _accessor(doc, blob, prim["attributes"]["NORMAL"])
                    # Normals transform by the inverse transpose.
                    n = n @ np.linalg.inv(world[:3, :3])
                    normals = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
                if "indices" in prim:
                    faces = _accessor(doc, blob, prim["indices"]).reshape(-1, 3)
                else:
                    faces = np.arange(len(pos)).reshape(-1, 3)
                mat = materials[prim["material"]] if "material" in prim else {}
                rgba = mat.get("pbrMetallicRoughness", {}).get(
                    "baseColorFactor", [0.6, 0.6, 0.6, 1.0])
                out.append(Primitive(material=mat.get("name", ""),
                                     colour=(float(rgba[0]), float(rgba[1]), float(rgba[2])),
                                     vertices=verts, faces=faces, normals=normals))
        for child in node.get("children", []):
            walk(child, world)

    scene = doc["scenes"][doc.get("scene", 0)]
    for root in scene["nodes"]:
        walk(root, np.eye(4))
    return out


__all__ = ["Primitive", "read_glb"]
