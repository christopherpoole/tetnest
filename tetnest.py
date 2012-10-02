import os
import sys
import math
import copy
import itertools
import logging

from collections import defaultdict
from pprint import pprint
from os import system

import numpy


class PLYMesh(object):
    def __init__(self, fn):
        self.filename = fn
        
        self.VERTEX_COUNT = 0
        self.FACE_COUNT = 0
        
        self.IS_PLY_FILE = False
        self.IS_ASCII = False
        
        self.header = []
        self.faces = []
        self.verts = []
        self.g4material = None
        
        self.ply_file = open(self.filename)
        
        for line in self.ply_file:
            line = line.strip()  
            self.header.append(line)
                        
            # If the first line of the file is not 'ply', we do not have a
            # Stanford polygon file format.
            if not self.IS_PLY_FILE:
                if line == 'ply':
                    self.IS_PLY_FILE = True
                    continue
                else:
                    raise IOError("%s is not a Stanford polygon file, the "\
                                  "first line must be 'ply'" % self.filename)
                    break
                    
            # If the second line of the file is not 'format ascii 1.0', we do
            # not have an ASCII file.      
            if not self.IS_ASCII:
                if line == 'format ascii 1.0':
                    self.IS_ASCII = True
                    continue
                else:
                    raise IOError("%s is not an ASCII file, the second line "\
                                  "must be 'format ascii 1.0'" % self.filename)
                    break
                       
            if line == 'end_header':
                break

            # Extract vertex and face counts
            parameter = line.split(' ')
            if parameter[0] == 'element':
                if parameter[1] == 'vertex':
                    self.VERTEX_COUNT = int(parameter[2])
                if parameter[1] == 'face':
                    self.FACE_COUNT = int(parameter[2])
           
            if parameter[0] == 'comment':
                if parameter[1] == 'g4material':
                    self.g4material = parameter[2]
        
        # If we parse the entire file without finding 'end_header' raise error.
        if self.header[-1] != 'end_header':
            raise IOError("%s does not have a defined header; 'end_header' "\
                          "missing" % self.filename)
           
        self.get_faces()
        self.get_verts()
        
    def get_faces(self):
        if len(self.faces) == 0:
            # The verts appear at the top of the file, so make sure
            # we get them first.
            self.get_verts()
            for i, line in enumerate(self.ply_file):
                #if i == self.FACE_COUNT:
                #    break
                line = line.strip().split(' ')
                # Ignore the first element, it is an vertex count. Let us assume
                # all faces are triangles only.
                self.faces.append(line[1:4])
        
        if len(self.faces) != self.FACE_COUNT:
            raise IOError("%s defines %i faces whilst %i are expected"
                           % (self.filename, len(self.faces), self.FACE_COUNT))
        
        self.faces = numpy.array(self.faces, dtype=float)
        return self.faces
        
        
    def get_verts(self):
        if len(self.verts) == 0:
            for i, line in enumerate(self.ply_file):
                line = line.strip().split(' ')
                self.verts.append(line)
                if i == self.VERTEX_COUNT-1:
                    break

        if len(self.verts) != self.VERTEX_COUNT:
            raise IOError("%s defines %i verts whilst %i are expected"
                           % (self.filename, len(self.faces), self.FACE_COUNT))
                           
        self.verts = numpy.array(self.verts, dtype=float)     
        return self.verts


class Mesh(object):
    def __init__(self):
        self.meshes = []
        self.vert_count = 0
        self.face_count = 0
        self.region_count = 0
        
        self.meshes = []
        self.regions = []
    
    def add(self, mesh, region):
        self.meshes.append(mesh)
        self.vert_count += len(mesh.verts)
        self.face_count += len(mesh.faces)
        
        self.regions.append(region)
        self.region_count += len(region)        
        
    def write(self, filename):
        f = file("%s.smesh" % filename, 'w')
        f.write("#Part 1 - nodes\n%i 3 0 0\n" % self.vert_count)
        
        counter = 0
        for m in self.meshes:
            for v in m.verts:
                f.write("%i %f %f %f\n" % (counter, v[0], v[1], v[2]))
                counter += 1
        assert counter == self.vert_count, "Vert count error: %i, %i" % (self.vert_count, counter)
        
        f.write("#Part 2 - faces\n%i 1\n" % self.face_count)
        counter = 0
        offset = 0
        for i, m in enumerate(self.meshes):
            if i > 0: offset += len(self.meshes[i-1].verts)
            for p in m.faces:
                f.write("3 %i %i %i %i\n" % (p[0]+offset, p[1]+offset, p[2]+offset, i))
                counter += 1
        assert counter == self.face_count, "Face count error: %i, %i" % (self.face_count, counter)
        
        f.write("#Part 3\n0\n")
        f.write("#Part 4\n%i\n" % self.region_count)
        counter = 0
        for i, r in enumerate(self.regions):
            for p in r:
                f.write("%i %f %f %f %i\n" % (counter, p[0], p[1], p[2], -i))
                counter += 1
                        
        f.close()
        

class TetMesh(object):
    """ Read a *.ele tetrahedron file and corresponding *.node file.
    """
    def __init__(self, name, read=True):
        self.file_name = name
        self.tets = []
        self.verts = []
        if read is True:
            self.read()            
            
    def read(self):           

        f = file("%s.ele" % self.file_name)
        self.tets = [map(int, l.split()[1:]) for l in f if l[0] != '#'][1:]
        

        f = file("%s.node" % self.file_name)
        self.verts = [map(float, l.split()[1:]) for l in f if l[0] != '#'][1:]
        self.verts_x = [v[0] for v in self.verts]
        self.verts_y = [v[1] for v in self.verts]
        self.verts_z = [v[2] for v in self.verts]
        

class MeshDiff(object):
    """ Search for the `target` Mesh in the `reference` Mesh and repair where
        the meshes to not correspond to the limit of `match_threshold`.
    """
    def __init__(self, target, reference, match_threshold=1, repair=True):
        self.target = target
        self.reference = reference
        self.match_threshold = match_threshold
               
        self.missing = self._check(self.target.verts, self.reference.verts)
        self.matches, self.match_distance = self._match(self.missing, self.reference.verts)
        self.valid_matches = self._check(self.matches, self.target.verts)



    def _check(self, target, reference):
        """ Check if any verts in `target` are missing from `reference`.
        """
        missing = []
        for v in target:
            try:
                reference.index(v)
            except:
                missing.append(tuple(v))
        return missing
    def _match(self, missing, reference):
        """ Search `reference` within `match_threshold` for the `missing`
            verts.
        """
        matches = []
        match_dist = []
        for m in missing:
            match = None
            match_threshold = self.match_threshold
            for t in reference:
                diff = map(lambda a, b: a - b, t, m)
                dist = math.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
                if dist < self.match_threshold and dist < match_threshold:
                    match = t
                    match_threshold = dist
            if match is not None:
                matches.append(tuple(match))
                match_dist.append(match_threshold)
        return matches, match_dist


class SplitMesh(object):
    def __init__(self, filename):
        self.filename = filename
        
        f = open("%s.ele" % self.filename).read().split('\n')  
        self.tets = defaultdict(list)
        for line in f[1:-2]:
            l = line.split()
            self.tets[l[5]].append(map(int, l[0:5]))

    def write(self):
        for k, v in self.tets.iteritems():
            f = file("%s_%s_.ele" % (self.filename, k), 'w')
            f.write("%i 4 0\n" % len(v))
            for t in enumerate(v):
                f.write("%s %s %s %s %s\n" % (t[0], t[1][1], t[1][2], t[1][3], t[1][4]))
            f.close()
            os.system("cp %s.node %s_%s_.node" % (self.filename, self.filename, k))


if __name__ == "__main__":
    names = sys.argv[1:]

    mesh = Mesh()
   
    # Nest Meshes 
    for n in names:
        ply_mesh = PLYMesh("%s" % n)
        os.system("tetgen -Y %s" % n)
        n = n[:-4] # assume .ply
        tet_mesh = TetMesh("%s.1" % n)
        
        verts = [tet_mesh.verts[tet_mesh.tets[0][v]] for v in range(0, 4)]
        x = min([v[0] for v in verts]) + (max([v[0] for v in verts]) \
                 - min([v[0] for v in verts]))/2
        y = min([v[1] for v in verts]) + (max([v[1] for v in verts]) \
                 - min([v[1] for v in verts]))/2
        z = min([v[2] for v in verts]) + (max([v[2] for v in verts]) \
                 - min([v[2] for v in verts]))/2
        mesh.add(ply_mesh, [(x, y, z)])
        
    mesh.write("combined")
    
    os.system("tetgen -YAO combined.smesh")
   
    # Split TetMeshes 
    split_mesh = SplitMesh("combined.1")
    split_mesh.write()
    
