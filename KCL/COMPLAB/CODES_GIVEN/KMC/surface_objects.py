# Version 1.0 by Mattia Fiorentini, January 2016
class Deposition(object):
    '''Method of deposition'''
    def __init__(self,drops,dep_rate,temperature):
        self.rate = dep_rate # sec E-16
        self.t = temperature # K
        self.drops = drops # units

class Substrate(object):
    '''Any surface on which droplet can diffuse'''
    def __init__(self,delta_surface,energy_ev_surface):
        self.de_diff = delta_surface # eV
        self.de_desorb = energy_ev_surface # eV
        
class Droplet(Substrate):
    '''Droplet diffusion properties'''
    def __init__(self,delta_stack,energy_ev_stack,energy_bind,delta_climb):
        super(Droplet, self).__init__(delta_stack,energy_ev_stack)
        self.de_bind = energy_bind # eV
        self.de_climb = delta_climb # eV
        
class Surface(Substrate):
    '''Base substrate of diffusion with periodic boundaries'''
    def __init__(self,sizes,delta_surface,energy_ev_surface):
        super(Surface, self).__init__(delta_surface,energy_ev_surface)
        self.lx = sizes[0] # unit lenght
        self.ly = sizes[1] # unit lenght
        self.max_lz = sizes[2] # unit lenght
        self.size = sizes
        self.adsorbate = {} # it is a dictionary: keys are _Coordiantes type objs and values are _Site type objs, it works as an hash table
        self.nads = 0
    
    def __iter__(self):
        #for v in self.adsorbate.values():
        #    yield v
        return iter(self.adsorbate.values()) # when an iterator-type object is needed, the list of values in the surface "dict" are given
    
    def add_part(self,s,part): # add a particle given a site obj
        if s._occ == self.max_lz: return None
        #print("Adding part. on site",s.pos)
        s._occ += 1
        s.de_diff = part.de_diff
        s.de_desorb = part.de_desorb
        self.adsorbate[s.coord] = s
    
    def remove_part(self,s): # remove a particle given a site obj
        #print("Removing part. on site",s.pos)
        if s._occ <= 1:
            self.adsorbate.pop(s.coord,None)
        else:
            s._occ -= 1
            self.adsorbate[s.coord] = s

    def remove_site(self,s):
        if s._occ > 0 :
            s._occ = 0
            self.adsorbate.pop(s.coord,None)
        
    def get_site(self,pos): # given a position, return the corresponding site obj
        coord = _Coordinate(self,pos)
        s = self.adsorbate.get(coord, _Site(self,pos)) # return site 
        return s
        
    def site_neighbours(self,s): # return neighbouring site-objects given a site
        neighbours = []
        for d in [0,1]:
            for v in [-1,1]:
                n_pos = s.pos[:]
                #print("n_pos1",n_pos)
                n_pos[d] += v
                #print("n_pos2",n_pos)
                neighbour = self.get_site(n_pos)
                neighbours.append(neighbour)
        return neighbours

class _Site(Substrate):
    '''State of a site belonging to a Surface object'''
    def __init__(self,surf,coords): 
        super(_Site,self).__init__(surf.de_diff,surf.de_desorb)
        self.coord = _Coordinate(surf,coords) 
        self.pos = self.coord.get_coord()
        self.de_diff = surf.de_diff
        self.de_desorb =surf.de_desorb
        self._occ = 0

    @property
    def occ(self):
        """Site occupation"""
        return self._occ

class _Coordinate(object):
    '''Position (given periodic boundary conditions) of a _Site object embedded on a Surface object, defined in term of non-periodic cartesian coordinates'''
    def __init__(self,surf,pos):
        self._bound = surf.size[0:2]
        # boundary condition according to std xy cart coord
        self.coord = pos[:]
        for i in range(2):
            while self.coord[i] < 0:
                self.coord[i] += self._bound[i]
            self.coord[i] = self.coord[i]%self._bound[i]
        self.x = self.coord[0]
        self.y = self.coord[1]

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
    
    def __hash__(self):
        return self.x + self.y*self._bound[0]
    
    def get_coord(self):
        #print("get_coord",self.coord)
        return self.coord    

class Move(object):
    '''Data structure to store the kinematics of particle movement'''
    def __init__(self,rate,in_surf,in_site,fi_surf,fi_site):
        self.rate = rate
        self.in_surf = in_surf
        self.in_site = in_site
        self.fi_surf = fi_surf
        self.fi_site = fi_site
