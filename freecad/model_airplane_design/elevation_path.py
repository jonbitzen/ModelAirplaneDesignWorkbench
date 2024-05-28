import FreeCAD as App
import numpy
import Part
from . import utilities
from typing import List

class Pose():
    def __init__(self, position: App.Vector, direction: App.Vector) -> None:
        self.position = position
        self.direction = direction

class PathInterval():
    def __init__(self, offset: float, edge: Part.Edge) -> None:
        self.edge = edge
        self.offset = offset
        # TODO: is there a universe where FirstParameter isn't zero, if so
        #       there is a discontinuity to cope with
        self.start = edge.FirstParameter + offset
        self.end = edge.LastParameter + offset

    def contains(self, location: float) -> bool:
        return location >= self.start and location <= self.end

    def get_pose_at(self, location: float) -> Pose:
        if not self.contains(location):
            return None
        loc: float = location - self.offset
        pos: App.Vector = self.edge.valueAt(loc)
        tan: App.Vector = self.edge.tangentAt(loc)
        return Pose(pos, tan)

class PathHelper():
    def __init__(self, path: List[Part.Edge]) -> None:

        self.path: List[PathInterval] = []
        self.length: float = 0

        origin: App.Vector = App.Vector(0,0,0)
        path.sort(key=lambda edge: edge.valueAt(edge.FirstParameter).distanceToPoint(origin))

        self.length = 0
        for edge in path:
            path_segment = \
                PathInterval(
                    self.length,
                    edge
                )

            self.path.append(path_segment)
            self.length += edge.Length

    def get_poses(self, num_poses: int) -> List[Pose]:
        pose_list: List[Pose] = []
        for edge_dist in numpy.linspace(0.0, self.length, num_poses):
            pose_generated: bool = False
            for edge_int in self.path:
                
                rib_pose: Pose = edge_int.get_pose_at(edge_dist)
                if rib_pose is not None:
                    pose_list.append(rib_pose)
                    pose_generated = True
                    continue

            if not pose_generated:
                print("PathHelper.get_poses: failed to generate pose at " +  str(edge_dist))
        return pose_list 