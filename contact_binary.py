import bpy
from bpy.props import FloatProperty
from math import *


G                  = 6.67384E-11  # Nm2kg-2
STELLAR_RADIUS     = 6.96342E8    # m
STELLAR_MASS       = 1.98855E30   # kg
STELLAR_LUMINOSITY = 3.846E26     # W 

def solve(fn, fn_dash, x=0, iteration_limit=20, accuracy=1):
    print('solve')
    y = fn(x)
    while iteration_limit > 0 and (y < -accuracy or y > accuracy):
        y_dash = fn_dash(x)
        print('- x=', x, ' y=', y, ' y\'=', y_dash)
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


class MeshGenerator:

    def __init__(self):
        self.points = []
        self.faces = []
            
    def add_point(self, x, y, z):
        n = len(self.points)
        p = [x, y, z]
        self.points.append(p)
        return n
    
    def add_row_points(self, x, n, offset):
        print('add_row_points(', x, ', ', n, ', ', offset)
        row = []
        
        max_z = solve(
            lambda z : self.surface(x, 0, z), \
            lambda z : self.dp_dz(x, 0, z), \
            self.radius)
            
        total = 5*n
        i = 0
        while i < total:
            angle = pi*(2*i+offset) / total
            z = max_z * cos(angle)
            y_guess = max_z * sin(angle)
            if y_guess < -1 or y_guess > 1:
                y = solve(
                    lambda y : self.surface(x, y, z), \
                    lambda y : self.dp_dy(x, y, z), \
                    y_guess)
            else:
                y = 0
            row.append(self.add_point(x / STELLAR_RADIUS, y / STELLAR_RADIUS, z / STELLAR_RADIUS))
            i += 1
       
        print()   
        print(row)
       
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
                  

class RotatingStar(MeshGenerator):

    def __init__(self, mass=STELLAR_MASS, radius=STELLAR_RADIUS, angular_speed=0):

        super(RotatingStar, self).__init__()
        
        if (angular_speed > 0):
            stationary_orbit = pow(G*mass/(angular_speed*angular_speed), 1.0/3.0)
            if (stationary_orbit <= 1.5 * radius):
                raise ValueError('Star will fly apart due to high rotation speed')
                
        self.mass = Mass(mass, 0, 0, 0)
        self.rotation = Rotation(angular_speed, 0, 0, 0)
        self.radius = radius
        
        self.surface_potential = self.potential(0, 0, radius)

    def potential(self, x, y, z):
        return self.mass.potential(x, y, z) + self.rotation.potential(x, y, z)

    def dp_dx(self, x, y, z):
        return self.mass.dp_dx(x, y, z) + self.rotation.dp_dx(x, y, z)

    def dp_dy(self, x, y, z):
        return self.mass.dp_dy(x, y, z) + self.rotation.dp_dy(x, y, z)

    def dp_dz(self, x, y, z):
        return self.mass.dp_dz(x, y, z) + self.rotation.dp_dz(x, y, z)

    def surface(self, x, y, z):
        return self.potential(x, y, z) - self.surface_potential

    def make_mesh(self, n):
        
        max_x = solve(
            lambda x : self.surface(x, 0, 0), \
            lambda x : self.dp_dx(x, 0, 0), \
            self.radius)

        row0 = [self.add_point(-max_x / STELLAR_RADIUS, 0, 0)]
        for i in range(1, n+1):
            x = -max_x * cos(i*pi / (3*n))
            row1 = self.add_row_points(x, i, 0)
            self.add_row_faces(row0, row1)
            row0 = row1
            
        for i in range(n+1, 2*n+1):
            x = -max_x * cos(i*pi / (3*n))
            row1 = self.add_row_points(x, n, i-n)
            self.add_row_faces(row0, row1)
            row0 = row1
            
        for i in range(2*n+1, 3*n):
            x = -max_x * cos(i*pi / (3*n))
            row1 = self.add_row_points(x, 3*n-i, 3*n-i)
            self.add_row_faces(row0, row1)
            row0 = row1
        row1 = [self.add_point(max_x / STELLAR_RADIUS, 0, 0)]
        self.add_row_faces(row0, row1) 
            
        return (self.points, self.faces)
    
    
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

#    rs = RotatingStar(angular_speed=0.0000025) # Sol
#    rs = RotatingStar(mass=6.7*STELLAR_MASS, radius=7.3*STELLAR_RADIUS, angular_speed=0.0000399) # Achernar
    rs = RotatingStar(mass=3.8*STELLAR_MASS, radius=2.5*STELLAR_RADIUS, angular_speed=0.000162) # Regulus
    (points, faces) = rs.make_mesh(10)
    add_mesh("Star", points, faces)
