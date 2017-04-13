import bpy
from bpy.props import FloatProperty
from math import *


G                  = 6.67384E-11  # Nm2kg-2
STELLAR_RADIUS     = 6.96342E8    # m
STELLAR_MASS       = 1.98855E30   # kg
STELLAR_LUMINOSITY = 3.846E26     # W 

MASS_LOSS_FACTOR   = 1.001

def solve(fn, fn_dash, x=0, iteration_limit=20, accuracy=1):
#    print('solve')
    y = fn(x)
    while iteration_limit > 0 and (y < -accuracy or y > accuracy):
        y_dash = fn_dash(x)
#        print('- x=', x, ' y=', y, ' y\'=', y_dash)
        x = x - y / y_dash
        y = fn(x)
        iteration_limit -= 1
    return x


class Mass:

    def __init__(self, mass = STELLAR_MASS, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass

    def potential(self, x, y, z):
        deltaX = x - self.x
        deltaY = y - self.y
        deltaZ = z - self.z
        return -G * self.mass / sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ)

    def dp_dx(self, x, y, z):
        deltaX = x - self.x
        deltaY = y - self.y
        deltaZ = z - self.z
        return G * self.mass * deltaX / pow(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ, 1.5)

    def dp_dy(self, x, y, z):
        deltaX = x - self.x
        deltaY = y - self.y
        deltaZ = z - self.z
        return G * self.mass * deltaY / pow(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ, 1.5)

    def dp_dz(self, x, y, z):
        deltaX = x - self.x
        deltaY = y - self.y
        deltaZ = z - self.z
        return G * self.mass * deltaZ / pow(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ, 1.5)
    
    def d2p_dx2(self, x, y, z):
        deltaX = x - self.x
        deltaY = y - self.y
        deltaZ = z - self.z
        return -2.0 * G * self.mass * (deltaX*deltaX - deltaY*deltaY - deltaZ*deltaZ) \
            / pow(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ, 2.5)


class Rotation:
    
    def __init__(self, angular_speed = 0.0, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z
        self.angular_speed = angular_speed

    def potential(self, x, y, z):
        deltaX = x - self.x
        deltaY = y - self.y
        return -0.5 * self.angular_speed * self.angular_speed * (deltaX*deltaX + deltaY*deltaY)

    def dp_dx(self, x, y, z):
        deltaX = x - self.x
        return -self.angular_speed * self.angular_speed * deltaX

    def dp_dy(self, x, y, z):
        deltaY = y - self.y
        return -self.angular_speed * self.angular_speed * deltaY

    def dp_dz(self, x, y, z):
        return 0
    
    def d2p_dx2(self, x, y, z):
        return -self.angular_speed * self.angular_speed
        

class MeshGenerator:

    def __init__(self):
        self.points = []
        self.faces = []
        self.offset = 0
            
    def add_point(self, x, y, z):
        n = len(self.points)
        p = [x, y, z]
        self.points.append(p)
        return n
    
    def add_row_points(self, x, rz, n, offset=0):
        row = []
        self.offset += offset
        
        max_z = solve(
            lambda z : self.potential(x, 0, z) - self.surface_potential, \
            lambda z : self.dp_dz(x, 0, z), \
            rz)
            
        total = 5*n
        i = 0
        while i < total:
            angle = self.offset + 2*pi*i / total
            z = max_z * cos(angle)
            y_guess = 0.99 * max_z * sin(angle)
            if y_guess < -1 or y_guess > 1:
                y = solve(
                    lambda y : self.potential(x, y, z) - self.surface_potential, \
                    lambda y : self.dp_dy(x, y, z), \
                    y_guess)
            else:
                y = 0
            row.append(self.add_point(x / STELLAR_RADIUS, y / STELLAR_RADIUS, z / STELLAR_RADIUS))
            i += 1
       
        return row

    def add_row_faces(self, row0, row1):
#        print('add_row_faces len(row0) = ',len(row0), ' len(row1) = ',len(row1))
        m = len(row0)
        n = len(row1)
        if m > 1: mc = m
        else: mc = 0
        if n > 1: nc = n
        else: nc = 0
        i = 0
        j = 0
        while i < mc or j < nc:
            ip = (2*i+1) * nc
            jp = (2*j+1) * mc
#            print('- i=', i, ' j=', j, ' ip=', ip, ' jp=', jp)
            if ip > jp:
                face = [row0[i%m], row1[j%n], row1[(j+1)%n]]
                j += 1
            else:
                face = [row0[i%m], row1[j%n], row0[(i+1)%m]]
                i += 1
 #           print('  - ', face)
            self.faces.append(face)
            

    def open_cone(self, x0, rx, rz, n):
        self.row0 = [self.add_point((x0-rx) / STELLAR_RADIUS, 0, 0)]
        for i in range(1, n+1):
            theta = i*pi / (3*n)
            x = x0 - rx * cos(theta)
            row1 = self.add_row_points(x, rz, i)
            self.add_row_faces(self.row0, row1)
            self.row0 = row1
            
    def cylinder(self, x0, x1, rz, n, offset=1):
        for i in range(1, n+1):
            x = ((n-i)*x0 + i*x1) / n
            row1 = self.add_row_points(x, rz, n, pi/(5*n))
            self.add_row_faces(self.row0, row1)
            self.row0 = row1
            
    def close_cone(self, x0, rx, rz, n):
        for i in range(1, n):
            theta = (2*n + i) * pi / (3*n)
            x = x0 - rx * cos(theta)
            row1 = self.add_row_points(x, rz, n-i)
            self.add_row_faces(self.row0, row1)
            self.row0 = row1
        row1 = [self.add_point((x0+rx) / STELLAR_RADIUS, 0, 0)]
        self.add_row_faces(self.row0, row1)
        self.row0 = None
                  

class RotatingStar(MeshGenerator):

    def __init__(self, mass=STELLAR_MASS, radius=STELLAR_RADIUS, angular_speed=0):

        super(RotatingStar, self).__init__()
                        
        self.mass = Mass(mass)
        self.rotation = Rotation(angular_speed)
        self.initial_estimate = radius
        
        self.surface_potential = self.potential(0, 0, radius)
        if angular_speed > 0:
            synchronous_orbit = pow(G*mass/(angular_speed*angular_speed), 1.0/3.0)
            synchronous_potential = self.potential(synchronous_orbit, 0, 0)
            if self.surface_potential > MASS_LOSS_FACTOR * synchronous_potential:
                print('Mass lost at equator due to high rotation speed')
                self.surface_potential = MASS_LOSS_FACTOR * synchronous_potential
                self.initial_estimate = synchronous_orbit * 0.6

    def potential(self, x, y, z):
        return self.mass.potential(x, y, z) + self.rotation.potential(x, y, z)

    def dp_dx(self, x, y, z):
        return self.mass.dp_dx(x, y, z) + self.rotation.dp_dx(x, y, z)

    def dp_dy(self, x, y, z):
        return self.mass.dp_dy(x, y, z) + self.rotation.dp_dy(x, y, z)

    def dp_dz(self, x, y, z):
        return self.mass.dp_dz(x, y, z) + self.rotation.dp_dz(x, y, z)

    def make_mesh(self, n):
        
        max_x = solve(
            lambda x : self.potential(x, 0, 0) - self.surface_potential, \
            lambda x : self.dp_dx(x, 0, 0), \
            self.initial_estimate)
            
        self.open_cone(0, max_x, self.initial_estimate, n)     
        self.cylinder(-max_x * cos(pi/3), max_x * cos(pi/3), self.initial_estimate, n)
        self.close_cone(0, max_x, self.initial_estimate, n)
            
        return (self.points, self.faces)
    
    
class BinaryStar(MeshGenerator):
    
    def __init__(self, \
            mass1=STELLAR_MASS, radius1=STELLAR_RADIUS, \
            mass2=STELLAR_MASS, radius2=STELLAR_RADIUS, \
            angular_speed=0.0001):

        super(BinaryStar, self).__init__()

        separation = pow(G*(mass1+mass2) / (angular_speed*angular_speed), 1.0/3.0)
        self.mass1 = Mass(mass1, x=-separation*mass2 / (mass1+mass2))
        self.initial_estimate1 = radius1
        self.mass2 = Mass(mass2, x=separation*mass1 / (mass1+mass2))
        self.initial_estimate2 = radius2
        self.rotation = Rotation(angular_speed, 0, 0, 0)
        
        self.surface_potential1 = self.potential(self.mass1.x, 0, radius1)
        self.surface_potential2 = self.potential(self.mass2.x, 0, radius2)

        self.__find_lagrange_points()
        self.__adjust_surface_potentials()
        
    def make_mesh(self, n):
        
        if self.common_envelope:
            self.surface_potential = self.surface_potential1
            self.__make_extended_spheroid(n)
            
        else:
            self.surface_potential = self.surface_potential1
            self.__make_spheroid(self.mass1.x, self.initial_estimate1, n)
            
            self.surface_potential = self.surface_potential2
            self.__make_spheroid(self.mass2.x, self.initial_estimate2, n)
                        
        return (self.points, self.faces)

    def __make_spheroid(self, x0, initial_estimate, n):
        min_x = x0 - solve(
            lambda x : self.potential(x, 0, 0) - self.surface_potential, \
            lambda x : self.dp_dx(x, 0, 0), \
            x0 - initial_estimate)
        max_x = solve(
            lambda x : self.potential(x, 0, 0) -self.surface_potential, \
            lambda x : self.dp_dx(x, 0, 0), \
            x0 + initial_estimate) - x0
            
        self.open_cone(x0, min_x, initial_estimate, n)     
        self.cylinder(x0-min_x * cos(pi/3), x0+max_x * cos(pi/3), initial_estimate, n)
        self.close_cone(x0, max_x, initial_estimate, n)
        
    def __make_extended_spheroid(self, n):
        min_x = self.mass1.x - solve(
            lambda x : self.potential(x, 0, 0) - self.surface_potential, \
            lambda x : self.dp_dx(x, 0, 0), \
            self.mass1.x - self.initial_estimate1)
        max_x = solve(
            lambda x : self.potential(x, 0, 0) - self.surface_potential, \
            lambda x : self.dp_dx(x, 0, 0), \
            self.mass2.x + self.initial_estimate2) - self.mass2.x

        m1_l1_x = (self.mass1.x + 2*self.l1x) / 3
        l1_m2_x = (2*self.l1x + self.mass2.x) / 3
        self.open_cone(self.mass1.x, min_x, self.initial_estimate1, n)     
        self.cylinder(self.mass1.x-min_x * cos(pi/3), m1_l1_x, self.initial_estimate1, n)
        self.cylinder(m1_l1_x, l1_m2_x, self.initial_estimate1, n)
        self.cylinder(l1_m2_x, self.mass2.x+max_x * cos(pi/3), self.initial_estimate1, n)
        self.close_cone(self.mass2.x, max_x, self.initial_estimate2, n)
    
    def __find_lagrange_points(self):
        self.l1x = solve(
            lambda x : self.dp_dx(x, 0, 0), \
            lambda x : self.d2p_dx2(x, 0, 0), \
            0)
        self.l1_potential = self.potential(self.l1x, 0, 0)

        self.l2x = solve(
            lambda x : self.dp_dx(x, 0, 0), \
            lambda x : self.d2p_dx2(x, 0, 0), \
            2.0*self.mass1.x - self.l1x)
        self.l2_potential = self.potential(self.l2x, 0, 0)

        self.l3x = solve(
            lambda x : self.dp_dx(x, 0, 0), \
            lambda x : self.d2p_dx2(x, 0, 0), \
            2.0*self.mass2.x - self.l1x)
        self.l3_potential = self.potential(self.l3x, 0, 0)
            
    def __adjust_surface_potentials(self):
        self.common_envelope = self.surface_potential1 > self.l1_potential and self.surface_potential2 > self.l1_potential
        if self.common_envelope:
            print('Stars share a common envelope')
            if self.surface_potential1 > self.surface_potential2:
                self.surface_potential1 = self.surface_potential2
            else:
                self.surface_potential2 = self.surface_potential1
                
        elif self.surface_potential1 > MASS_LOSS_FACTOR * self.l1_potential:
            print('Star A loses mass through L1 point')
            self.surface_potential1 = MASS_LOSS_FACTOR * self.l1_potential
                
        elif self.surface_potential2 > MASS_LOSS_FACTOR * self.l1_potential:
            print('Star B loses mass through L1 point')
            self.surface_potential2 = MASS_LOSS_FACTOR * self.l1_potential
            
        if self.common_envelope:
            if self.surface_potential1 > MASS_LOSS_FACTOR * self.l2_potential:
                print('Common envelope loses mass through L2 point')
            if self.surface_potential2 > MASS_LOSS_FACTOR * self.l3_potential:
                print('Common envelope loses mass through L3 point')
            self.surface_potential1 = min(self.surface_potential1, MASS_LOSS_FACTOR * self.l2_potential, MASS_LOSS_FACTOR * self.l3_potential)
            self.surface_potential2 = min(self.surface_potential2, MASS_LOSS_FACTOR * self.l2_potential, MASS_LOSS_FACTOR * self.l3_potential)

        while self.potential(self.mass1.x - self.initial_estimate1, 0, 0) > self.surface_potential1:
            print('reducing initial_estimate1')
            self.initial_estimate1 /= 2.0        
        while self.potential(self.mass2.x + self.initial_estimate2, 0, 0) > self.surface_potential2:
            print('reducing initial_estimate2')
            self.initial_estimate2 /= 2.0        

    def __str__(self):
        str = 'Star A: mass {:.2f} position {:.2f} surface {:f}\n'.format(
            self.mass1.mass / STELLAR_MASS, \
            self.mass1.x / STELLAR_RADIUS, \
            self.surface_potential1 / 1e9)
        str += 'Star B: mass {:.2f} position {:.2f} surface {:f}\n'.format(
            self.mass2.mass / STELLAR_MASS, \
            self.mass2.x / STELLAR_RADIUS, \
            self.surface_potential2 / 1e9)
        str += 'L1 position: {:.2f} potential {:f}\n'.format(
            self.l1x / STELLAR_RADIUS, \
            self.l1_potential / 1e9)
        str += 'L2 position: {:.2f} potential {:f}\n'.format(
            self.l2x / STELLAR_RADIUS, \
            self.l2_potential / 1e9)
        str += 'L3 position: {:.2f} potential {:f}\n'.format(
            self.l3x / STELLAR_RADIUS, \
            self.l3_potential / 1e9)
        return str

    def potential(self, x, y, z):
        return self.mass1.potential(x, y, z) + self.mass2.potential(x, y, z) + self.rotation.potential(x, y, z)

    def dp_dx(self, x, y, z):
        return self.mass1.dp_dx(x, y, z) + self.mass2.dp_dx(x, y, z) + self.rotation.dp_dx(x, y, z)

    def dp_dy(self, x, y, z):
        return self.mass1.dp_dy(x, y, z) + self.mass2.dp_dy(x, y, z) + self.rotation.dp_dy(x, y, z)

    def dp_dz(self, x, y, z):
        return self.mass1.dp_dz(x, y, z) + self.mass2.dp_dz(x, y, z) + self.rotation.dp_dz(x, y, z)

    def d2p_dx2(self, x, y, z):
        return self.mass1.d2p_dx2(x, y, z) + self.mass2.d2p_dx2(x, y, z) + self.rotation.d2p_dx2(x, y, z)
    
    
def add_mesh(name, points, faces):

    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(points, [], faces)
    mesh.update()

    obj = bpy.data.objects.new(name, mesh)
    obj.location = bpy.context.scene.cursor_location
    bpy.context.scene.objects.link(obj)


class StarOperator(bpy.types.Operator):
    bl_idname = "mesh.star_add"
    bl_label = "Star"
        
    mass = FloatProperty(
        name = "Mass",
        description = "Mass in solar masses",
        min = 0.0001,
        max = 1000,
        default = 1)
        
    radius = FloatProperty(
        name = "Radius",
        description = "Radius in solar radii",
        min = 0.0001,
        max = 1000,
        default = 1)
        
    spin = FloatProperty(
        name = "Spin",
        description = "Rotation period in hours (0 for no rotation)",
        min = 0,
        max = 1000,
        default = 24,
        precision = 3,
        step = 0.0001)
        
    def execute(self, context):
        angular_speed = spin
        if angular_speed != 0:
            angular_speed = 2 * pi * 3600 / angular_speed
        star = RotatingStar(mass * STELLAR_MASS, angular_speed)
        self.__make_row_points([], star.surface, 0, 1, 0)
#        add_mesh(None, None)
        return {'FINISHED'}
           
       
class BinaryStarOperator(bpy.types.Operator):
    bl_idname = "mesh.binary_star_add"
    bl_label = "Binary star"

    mass1 = FloatProperty(
        name = "Mass 1",
        min = 0.0001,
        max = 1000,
        default = 1)
        
    radius1 = FloatProperty(
        name = "Radius 1",
        min = 0.0001,
        max = 1000,
        default = 1)
        
    mass2 = FloatProperty(
        name = "Mass 2",
        min = 0.0001,
        max = 1000,
        default = 1)
        
    radius2 = FloatProperty(
        name = "Radius 2",
        min = 0.0001,
        max = 1000,
        default = 1)
        
    separation = FloatProperty(
        name = "Separation",
        min = 0.0001,
        max = 1000,
        default = 1)
        
    def execute(self, context):
        None
                    
        
class Stars_add_menu(bpy.types.Menu):
    bl_idname = "Stars_add_menu"
    bl_label = "Stars"
 
    def draw(self, context):
        self.layout.operator_context = 'INVOKE_REGION_WIN'
        self.layout.operator(StarOperator.bl_idname, text = "Star")
        self.layout.operator(BinaryStarOperator.bl_idname, text = "Binary star")
      
                    
def menu_func(self, context):
    self.layout.menu(Stars_add_menu.bl_idname, icon="PLUGIN")


def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_mesh_add.append(menu_func)


def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_mesh_add.remove(menu_func)


if __name__ == "__main__":
#    register()

#    s = RotatingStar(angular_speed=0.0000025) # Sol
#    s = RotatingStar(mass=6.7*STELLAR_MASS, radius=7.3*STELLAR_RADIUS, angular_speed=0.0000399) # Achernar
#    s = RotatingStar(mass=3.8*STELLAR_MASS, radius=2.5*STELLAR_RADIUS, angular_speed=0.000162) # Regulus
    s = BinaryStar(\
        mass1=0.95*STELLAR_MASS, radius1=0.88*STELLAR_RADIUS, \
        mass2=0.55*STELLAR_MASS, radius2=0.66*STELLAR_RADIUS, \
        angular_speed=0.00033) #angular_speed=0.0002716)
    (points, faces) = s.make_mesh(15)
    add_mesh("Star", points, faces)
    print(s)
