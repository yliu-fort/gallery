#include "isosurface.h"
#include <iostream>
#include <cmath>
#include <vector>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

//GLFW
#include <GLFW/glfw3.h>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#include "camera.h"
#include "shader.h"
#include "cmake_source_dir.h"

static int edgeTable[256] =
{
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

// triangle table maps same cube vertex index to a list of up to 5 triangles
// which are built from the interpolated edge vertices
#define X 255
static int triTable[256][16] =
{
    {X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {0, 8, 3, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {0, 1, 9, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {1, 8, 3, 9, 8, 1, X, X, X, X, X, X, X, X, X, X},
    {1, 2, 10, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {0, 8, 3, 1, 2, 10, X, X, X, X, X, X, X, X, X, X},
    {9, 2, 10, 0, 2, 9, X, X, X, X, X, X, X, X, X, X},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, X, X, X, X, X, X, X},
    {3, 11, 2, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {0, 11, 2, 8, 11, 0, X, X, X, X, X, X, X, X, X, X},
    {1, 9, 0, 2, 3, 11, X, X, X, X, X, X, X, X, X, X},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, X, X, X, X, X, X, X},
    {3, 10, 1, 11, 10, 3, X, X, X, X, X, X, X, X, X, X},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, X, X, X, X, X, X, X},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, X, X, X, X, X, X, X},
    {9, 8, 10, 10, 8, 11, X, X, X, X, X, X, X, X, X, X},
    {4, 7, 8, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {4, 3, 0, 7, 3, 4, X, X, X, X, X, X, X, X, X, X},
    {0, 1, 9, 8, 4, 7, X, X, X, X, X, X, X, X, X, X},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, X, X, X, X, X, X, X},
    {1, 2, 10, 8, 4, 7, X, X, X, X, X, X, X, X, X, X},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, X, X, X, X, X, X, X},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, X, X, X, X, X, X, X},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, X, X, X, X},
    {8, 4, 7, 3, 11, 2, X, X, X, X, X, X, X, X, X, X},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, X, X, X, X, X, X, X},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, X, X, X, X, X, X, X},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, X, X, X, X},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, X, X, X, X, X, X, X},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, X, X, X, X},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, X, X, X, X},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, X, X, X, X, X, X, X},
    {9, 5, 4, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {9, 5, 4, 0, 8, 3, X, X, X, X, X, X, X, X, X, X},
    {0, 5, 4, 1, 5, 0, X, X, X, X, X, X, X, X, X, X},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, X, X, X, X, X, X, X},
    {1, 2, 10, 9, 5, 4, X, X, X, X, X, X, X, X, X, X},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, X, X, X, X, X, X, X},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, X, X, X, X, X, X, X},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, X, X, X, X},
    {9, 5, 4, 2, 3, 11, X, X, X, X, X, X, X, X, X, X},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, X, X, X, X, X, X, X},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, X, X, X, X, X, X, X},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, X, X, X, X},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, X, X, X, X, X, X, X},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, X, X, X, X},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, X, X, X, X},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, X, X, X, X, X, X, X},
    {9, 7, 8, 5, 7, 9, X, X, X, X, X, X, X, X, X, X},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, X, X, X, X, X, X, X},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, X, X, X, X, X, X, X},
    {1, 5, 3, 3, 5, 7, X, X, X, X, X, X, X, X, X, X},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, X, X, X, X, X, X, X},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, X, X, X, X},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, X, X, X, X},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, X, X, X, X, X, X, X},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, X, X, X, X, X, X, X},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, X, X, X, X},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, X, X, X, X},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, X, X, X, X, X, X, X},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, X, X, X, X},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, X},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, X},
    {11, 10, 5, 7, 11, 5, X, X, X, X, X, X, X, X, X, X},
    {10, 6, 5, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {0, 8, 3, 5, 10, 6, X, X, X, X, X, X, X, X, X, X},
    {9, 0, 1, 5, 10, 6, X, X, X, X, X, X, X, X, X, X},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, X, X, X, X, X, X, X},
    {1, 6, 5, 2, 6, 1, X, X, X, X, X, X, X, X, X, X},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, X, X, X, X, X, X, X},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, X, X, X, X, X, X, X},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, X, X, X, X},
    {2, 3, 11, 10, 6, 5, X, X, X, X, X, X, X, X, X, X},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, X, X, X, X, X, X, X},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, X, X, X, X, X, X, X},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, X, X, X, X},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, X, X, X, X, X, X, X},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, X, X, X, X},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, X, X, X, X},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, X, X, X, X, X, X, X},
    {5, 10, 6, 4, 7, 8, X, X, X, X, X, X, X, X, X, X},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, X, X, X, X, X, X, X},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, X, X, X, X, X, X, X},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, X, X, X, X},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, X, X, X, X, X, X, X},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, X, X, X, X},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, X, X, X, X},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, X},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, X, X, X, X, X, X, X},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, X, X, X, X},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, X, X, X, X},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, X},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, X, X, X, X},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, X},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, X},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, X, X, X, X},
    {10, 4, 9, 6, 4, 10, X, X, X, X, X, X, X, X, X, X},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, X, X, X, X, X, X, X},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, X, X, X, X, X, X, X},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, X, X, X, X},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, X, X, X, X, X, X, X},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, X, X, X, X},
    {0, 2, 4, 4, 2, 6, X, X, X, X, X, X, X, X, X, X},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, X, X, X, X, X, X, X},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, X, X, X, X, X, X, X},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, X, X, X, X},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, X, X, X, X},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, X},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, X, X, X, X},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, X},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, X, X, X, X, X, X, X},
    {6, 4, 8, 11, 6, 8, X, X, X, X, X, X, X, X, X, X},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, X, X, X, X, X, X, X},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, X, X, X, X},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, X, X, X, X},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, X, X, X, X, X, X, X},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, X, X, X, X},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, X},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, X, X, X, X, X, X, X},
    {7, 3, 2, 6, 7, 2, X, X, X, X, X, X, X, X, X, X},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, X, X, X, X},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, X},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, X},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, X, X, X, X},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, X},
    {0, 9, 1, 11, 6, 7, X, X, X, X, X, X, X, X, X, X},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, X, X, X, X},
    {7, 11, 6, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {7, 6, 11, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {3, 0, 8, 11, 7, 6, X, X, X, X, X, X, X, X, X, X},
    {0, 1, 9, 11, 7, 6, X, X, X, X, X, X, X, X, X, X},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, X, X, X, X, X, X, X},
    {10, 1, 2, 6, 11, 7, X, X, X, X, X, X, X, X, X, X},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, X, X, X, X, X, X, X},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, X, X, X, X, X, X, X},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, X, X, X, X},
    {7, 2, 3, 6, 2, 7, X, X, X, X, X, X, X, X, X, X},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, X, X, X, X, X, X, X},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, X, X, X, X, X, X, X},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, X, X, X, X},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, X, X, X, X, X, X, X},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, X, X, X, X},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, X, X, X, X},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, X, X, X, X, X, X, X},
    {6, 8, 4, 11, 8, 6, X, X, X, X, X, X, X, X, X, X},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, X, X, X, X, X, X, X},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, X, X, X, X, X, X, X},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, X, X, X, X},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, X, X, X, X, X, X, X},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, X, X, X, X},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, X, X, X, X},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, X},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, X, X, X, X, X, X, X},
    {0, 4, 2, 4, 6, 2, X, X, X, X, X, X, X, X, X, X},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, X, X, X, X},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, X, X, X, X, X, X, X},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, X, X, X, X},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, X, X, X, X, X, X, X},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, X},
    {10, 9, 4, 6, 10, 4, X, X, X, X, X, X, X, X, X, X},
    {4, 9, 5, 7, 6, 11, X, X, X, X, X, X, X, X, X, X},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, X, X, X, X, X, X, X},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, X, X, X, X, X, X, X},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, X, X, X, X},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, X, X, X, X, X, X, X},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, X, X, X, X},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, X, X, X, X},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, X},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, X, X, X, X, X, X, X},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, X, X, X, X},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, X, X, X, X},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, X},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, X, X, X, X},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, X},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, X},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, X, X, X, X},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, X, X, X, X, X, X, X},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, X, X, X, X},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, X, X, X, X},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, X, X, X, X, X, X, X},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, X, X, X, X},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, X},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, X},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, X, X, X, X},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, X, X, X, X},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, X, X, X, X, X, X, X},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, X},
    {1, 5, 6, 2, 1, 6, X, X, X, X, X, X, X, X, X, X},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, X},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, X, X, X, X},
    {0, 3, 8, 5, 6, 10, X, X, X, X, X, X, X, X, X, X},
    {10, 5, 6, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {11, 5, 10, 7, 5, 11, X, X, X, X, X, X, X, X, X, X},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, X, X, X, X, X, X, X},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, X, X, X, X, X, X, X},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, X, X, X, X},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, X, X, X, X, X, X, X},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, X, X, X, X},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, X, X, X, X},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, X},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, X, X, X, X, X, X, X},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, X, X, X, X},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, X, X, X, X},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, X},
    {1, 3, 5, 3, 7, 5, X, X, X, X, X, X, X, X, X, X},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, X, X, X, X, X, X, X},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, X, X, X, X, X, X, X},
    {9, 8, 7, 5, 9, 7, X, X, X, X, X, X, X, X, X, X},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, X, X, X, X, X, X, X},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, X, X, X, X},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, X, X, X, X},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, X},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, X, X, X, X},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, X},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, X},
    {9, 4, 5, 2, 11, 3, X, X, X, X, X, X, X, X, X, X},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, X, X, X, X},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, X, X, X, X, X, X, X},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, X},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, X, X, X, X},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, X, X, X, X, X, X, X},
    {0, 4, 5, 1, 0, 5, X, X, X, X, X, X, X, X, X, X},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, X, X, X, X},
    {9, 4, 5, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, X, X, X, X, X, X, X},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, X, X, X, X},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, X, X, X, X},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, X},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, X, X, X, X},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, X},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, X, X, X, X, X, X, X},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, X, X, X, X},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, X, X, X, X},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, X},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, X},
    {1, 10, 2, 8, 7, 4, X, X, X, X, X, X, X, X, X, X},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, X, X, X, X, X, X, X},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, X, X, X, X},
    {4, 0, 3, 7, 4, 3, X, X, X, X, X, X, X, X, X, X},
    {4, 8, 7, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {9, 10, 8, 10, 11, 8, X, X, X, X, X, X, X, X, X, X},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, X, X, X, X, X, X, X},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, X, X, X, X, X, X, X},
    {3, 1, 10, 11, 3, 10, X, X, X, X, X, X, X, X, X, X},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, X, X, X, X, X, X, X},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, X, X, X, X},
    {0, 2, 11, 8, 0, 11, X, X, X, X, X, X, X, X, X, X},
    {3, 2, 11, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, X, X, X, X, X, X, X},
    {9, 10, 2, 0, 9, 2, X, X, X, X, X, X, X, X, X, X},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, X, X, X, X},
    {1, 10, 2, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {1, 3, 8, 9, 1, 8, X, X, X, X, X, X, X, X, X, X},
    {0, 9, 1, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {0, 3, 8, X, X, X, X, X, X, X, X, X, X, X, X, X},
    {X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X}
};
#undef X

// number of vertices for each case above
static int numVertsTable[256] =
{
    0,
    3,
    3,
    6,
    3,
    6,
    6,
    9,
    3,
    6,
    6,
    9,
    6,
    9,
    9,
    6,
    3,
    6,
    6,
    9,
    6,
    9,
    9,
    12,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    9,
    3,
    6,
    6,
    9,
    6,
    9,
    9,
    12,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    9,
    6,
    9,
    9,
    6,
    9,
    12,
    12,
    9,
    9,
    12,
    12,
    9,
    12,
    15,
    15,
    6,
    3,
    6,
    6,
    9,
    6,
    9,
    9,
    12,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    9,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    15,
    9,
    12,
    12,
    15,
    12,
    15,
    15,
    12,
    6,
    9,
    9,
    12,
    9,
    12,
    6,
    9,
    9,
    12,
    12,
    15,
    12,
    15,
    9,
    6,
    9,
    12,
    12,
    9,
    12,
    15,
    9,
    6,
    12,
    15,
    15,
    12,
    15,
    6,
    12,
    3,
    3,
    6,
    6,
    9,
    6,
    9,
    9,
    12,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    9,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    15,
    9,
    6,
    12,
    9,
    12,
    9,
    15,
    6,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    15,
    9,
    12,
    12,
    15,
    12,
    15,
    15,
    12,
    9,
    12,
    12,
    9,
    12,
    15,
    15,
    12,
    12,
    9,
    15,
    6,
    15,
    12,
    6,
    3,
    6,
    9,
    9,
    12,
    9,
    12,
    12,
    15,
    9,
    12,
    12,
    15,
    6,
    9,
    9,
    6,
    9,
    12,
    12,
    15,
    12,
    15,
    15,
    6,
    12,
    9,
    15,
    12,
    9,
    6,
    12,
    3,
    9,
    12,
    12,
    15,
    12,
    15,
    9,
    12,
    12,
    15,
    15,
    6,
    9,
    12,
    6,
    3,
    6,
    9,
    9,
    6,
    9,
    12,
    6,
    3,
    9,
    6,
    12,
    3,
    6,
    3,
    3,
    0,
};

// ------------------------------------------------------------------------
static unsigned int triTableTex;
static unsigned int numVertsTableTex;
static unsigned int dummyVAO;
static Shader _shaderHandle;
#define FACE_CULLING

void IsoSurface::ReloadShader()
{
    _shaderHandle.reload_shader_program_from_files(
                FP("glsl/mc.vert"),FP("glsl/mc.frag"),FP("glsl/mc.geom.glsl"));
}

void IsoSurface::Init()
{
    ReloadShader();

    //Triangle Table texture//
    //This texture store the vertex index list forgridPos
    //generating the triangles of each configurations.
    glGenTextures(1, &triTableTex);
    glActiveTexture(GL_TEXTURE10);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, triTableTex);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage1D( GL_TEXTURE_1D, 0, GL_R16I, 16*256, 0,
                  GL_RED_INTEGER, GL_INT, &triTable);

    //Numverts Table
    glGenTextures(1, &numVertsTableTex);
    glActiveTexture(GL_TEXTURE11);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, numVertsTableTex);
    //Integer textures must use nearest filtering mode
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    //We create an integer texture with new GL_EXT_texture_integer formats
    glTexImage1D( GL_TEXTURE_1D, 0, GL_R16I, 256, 0,
                  GL_RED_INTEGER, GL_INT, &numVertsTable);

    // Dummy VAO
    glGenVertexArrays(1, &dummyVAO);
}

template<typename T, typename S>
void IsoSurface::Draw(int target, float isoValue,
                      const T& camera,
                      const S& gridSize)
{

    // Set uniform attributes
    _shaderHandle.use();
    _shaderHandle.setMat4("projectionMatrix", camera.GetFrustumMatrix()*glm::scale(glm::mat4(1.0),
                                              glm::vec3(gridSize.x/(float)gridSize.y,1.0f,gridSize.z/(float)gridSize.y)));
    _shaderHandle.setVec3("viewPos",camera.Position);
    _shaderHandle.setInt("volumeTex",target);
    _shaderHandle.setInt("triTex",10);
    _shaderHandle.setInt("numVertsTex",11);
    _shaderHandle.setFloat("isoValue",isoValue);
    _shaderHandle.setVec3i("gridSize",gridSize);
    _shaderHandle.setVec3("voxelSize",1.0f/glm::vec3(gridSize));

    // Bind texture
    glActiveTexture(GL_TEXTURE10);
    glBindTexture(GL_TEXTURE_1D, triTableTex);
    glActiveTexture(GL_TEXTURE11);
    glBindTexture(GL_TEXTURE_1D, numVertsTableTex);

    // Draw
#ifdef FACE_CULLING
    glEnable(GL_CULL_FACE); // Enable to decrease computing load
    glCullFace(GL_FRONT);
#endif
    glBindVertexArray(dummyVAO);
    glDrawArraysInstanced(GL_POINTS,0,1,gridSize.x*gridSize.y*gridSize.z);
    glBindVertexArray(0);
#ifdef FACE_CULLING
    glDisable(GL_CULL_FACE); // Enable to decrease computing load
#endif
}

template void IsoSurface::Draw<Camera, glm::ivec3>(int, float, const Camera&, const glm::ivec3&);
//template void IsoSurface::Draw<double>(int, double, const glm::mat4&, const glm::ivec3&);
//template void IsoSurface::Draw<int>(int, int, const glm::mat4&, const glm::ivec3&);

// Initialize 3d datafield
//Datafield//
static unsigned int _demoVolumeTex;

// an interesting field function
static float torus(float x, float y, float z)
{
    float R2 = 0.4,a2 = 0.05;
    return (x*x + y*y + z*z + R2 - a2)*(x*x + y*y + z*z + R2 - a2) - 4.0*R2*(x*x + z*z)+1.0;
}
static float tangle(float x, float y, float z)
{
    x *= 3.0f;
    y *= 3.0f;
    z *= 3.0f;
    return (x*x*x*x - 5.0f*x*x +y*y*y*y - 5.0f*y*y +z*z*z*z - 5.0f*z*z + 11.8f) * 0.2f + 0.5f;
}
static float surf_of_evolution(float x, float z, float y)
{
    x *= 3.0f;
    y *= 3.0f;
    z *= 3.0f;
    return (x*x + y*y - (log(z+3.2)*log(z+3.2))-0.02 + 1.0);
}

void IsoSurface::Demo(int nx, int ny, int nz)
{
    //Store the volume data to polygonise
    glGenTextures(1, &_demoVolumeTex);
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, _demoVolumeTex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    //Generate a distance field to the center of the cube
    glm::vec4 *dataField=new glm::vec4[nx*ny*nz];
    for(int k=0; k<nz; k++)
        for(int j=0; j<ny; j++)
            for(int i=0; i<nx; i++){
                //dataField[i+j*nx+k*nx*ny].x = rand() / double(RAND_MAX);
                //dataField[i+j*nx+k*nx*ny].y = rand() / double(RAND_MAX);
                //dataField[i+j*nx+k*nx*ny].z = rand() / double(RAND_MAX);
                dataField[i+j*nx+k*nx*ny].x=glm::distance(glm::vec3(i, j, k),glm::vec3(64))/64.0f;
                dataField[i+j*nx+k*nx*ny].y=glm::distance(glm::vec3(i, j, k),glm::vec3(64))/64.0f;
                dataField[i+j*nx+k*nx*ny].z=glm::distance(glm::vec3(i, j, k),glm::vec3(64))/64.0f;
                //dataField[i+j*nx+k*nx*ny].w = sin(0.4*i)+sin(0.4*j) + sin(0.4*k);
                //dataField[i+j*nx+k*nx*ny].w = sin(0.4*i + 0.4*j + 0.4*k);
                //dataField[i+j*nx+k*nx*ny].w = sin(0.1*i * 0.1*j * 0.1*k);
                //dataField[i+j*nx+k*nx*ny].w = torus(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                dataField[i+j*nx+k*nx*ny].w = tangle(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                //dataField[i+j*nx+k*nx*ny].w = surf_of_evolution(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                //std::cout << dataField[i+j*nx+k*nx*ny] << std::endl;
            }
    glTexImage3D( GL_TEXTURE_3D, 0, GL_RGBA32F, nx, ny, nz, 0,
                  GL_RGBA, GL_FLOAT, dataField);
    delete [] dataField;
    dataField=NULL;
}
