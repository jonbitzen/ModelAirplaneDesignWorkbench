import FreeCAD as App
import numpy
import Part
from . import utilities
from typing import List

class Pose():
    """
        Wrapper for the position and tangent direction of a single pose generated
        along the path

        Parameters
        ----------
        
        position: App.Vector
            Cartesian location of the pose along the path

        direction: App.Vector
            Tangent vector of the pose along the path
    """
    def __init__(self, position: App.Vector, direction: App.Vector) -> None:
        self.position = position
        self.direction = direction

class PathInterval():
    """
        Encapsulates a single edge that is part of a multi-edge wing elevation 
        path, along with an offset that represents the starting point of this
        particular edge along the entire path.  This allows to calculate a
        tangent pose long the elevation path in absolute coordinates given a 
        path-length location

        Parameters
        ----------
        
        offset: float
            total distance from the beginning of the elevation path to the start
            point of the edge encapsulated by the PathInterval.

        edge: Part.Edge
            the edge that is encapsulated by the PathInterval
    """
    def __init__(self, offset: float, edge: Part.Edge) -> None:
        self.edge = edge
        self.offset = offset
        
        self.start = offset
        self.end = edge.Length + offset

        # alias the first and last parameter to find which is nearest the origin
        first_param_pt: App.Vector = edge.valueAt(edge.FirstParameter)
        second_param_pt: App.Vector = edge.valueAt(edge.LastParameter)

        first_pt_dist = first_param_pt.distanceToPoint(utilities.origin)
        second_pt_dist = second_param_pt.distanceToPoint(utilities.origin)

        self.direction_reversed = False
        if first_pt_dist > second_pt_dist:
            self.direction_reversed = True        

    def contains(self, location: float) -> bool:
        """
            Determines whether an absolute path-length location is contained
            within the current PathInterval

            Parameters
            ----------
            
            location: float
                The path-length location along the entire path to query

            bool:
                True if the desired path location is contained within this
                interval, False otherwise
        """
        return location >= self.start and location <= self.end

    def get_pose_at(self, location: float) -> Pose:
        """
            Calculates a tangent line pose located at a designated path-length
            location.

            Parameters
            ----------
            
            location: float
                a location along the total length of the elevation path at which
                to query for a tangent pose

            Returns
            -------

            Pose
                cartesian location along the PathInterval as well as a tangent 
                vector at that location, corresponding to the desired path-length
                location

                None if the desired path-length location is not in the current
                interval
        """
        if not self.contains(location):
            return None
        loc: float = location - self.offset

        if self.direction_reversed:
            loc = self.edge.Length - loc

        param_at_loc = self.edge.getParameterByLength(loc)

        pos: App.Vector = self.edge.valueAt(param_at_loc)
        tan: App.Vector = self.edge.tangentAt(param_at_loc)
 
        return Pose(pos, tan)

class PathHelper():
    """
        Wraps a list of edges that represent a wing elevation path.  The helper
        class allows you to locate tangent poses along N evenly spaced intervals
        along the path

        Parameters
        ----------
        path: List[Part.Edge]
            The list of edges that represent the path.  The edges should form a 
            continuous path, but may be in any order.  Endpoints for edges may 
            be either left-to-right or right-to-left, depending on how they were
            originally modelled
    """
    def __init__(self, path: List[Part.Edge]) -> None:

        self.path: List[PathInterval] = []
        self.length: float = 0

        # we want the edges to be in spatial order so they can be traversed from 
        # the wing root to the wing tip.
        path.sort(
            key=lambda edge: 
            min(
                edge.valueAt(edge.FirstParameter).distanceToPoint(utilities.origin),
                edge.valueAt(edge.LastParameter).distanceToPoint(utilities.origin)
            )
        )

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
        """
            Calculates a set of N evenly-spaced tangent poses along the entire
            length of the path defined by the set of edges

            Parameters
            ----------

            num_poses: int
                The number of poses to generate along the path

            Returns
            -------

            List[Pose]
                The list of generated poses along the path
        """
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
    
    def show(self, num_poses: int) -> None:
        pose_list = self.get_poses(num_poses)
        for pose in pose_list:
            utilities.draw_point(pose.position)
            utilities.draw_line(pose.position, pose.direction, line_len=10.0)