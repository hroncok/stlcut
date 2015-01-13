/* Copyright 2015 Miro Hrončok <miro@hroncok.cz>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */
#include <iostream>
#include <deque>
#include <set>
#include <math.h>
#include <admesh/stl.h>
#define TOLERANCE 0.000001

enum stl_position { above, on, below };

stl_vertex normalize(stl_vertex in) {
  double size = sqrt((double)in.x*in.x+(double)in.y*in.y+(double)in.z*in.z);
  in.x = in.x/size;
  in.y = in.y/size;
  in.z = in.z/size;
  return in;
}

float scalar(stl_vertex a, stl_vertex b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

struct stl_plane {
  float x;
  float y;
  float z;
  float d;
  stl_vertex a;
  stl_vertex b;
  
  stl_plane(float x, float y, float z, float d) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->d = d;
    
    if (x == 0 && y == 0) {
      a.x = 1; a.y = 0; a.z = 0;
      b.x = 0; b.y = 1; b.z = 0;
    } else if (y == 0 && z == 0) {
      a.x = 0; a.y = 1; a.z = 0;
      b.x = 0; b.y = 0; b.z = 1;
    } else if (x == 0 && z == 0) {
      a.x = 1; a.y = 0; a.z = 0;
      b.x = 0; b.y = 0; b.z = 1;
    } else {
      a.x = y; a.y = -x; a.z = 0;
      a = normalize(a);
      b.x = 0; b.y = z; b.z = -y;
      
      float r = scalar(a,b);
      b.x -= a.x*r;
      b.y -= a.y*r;
      b.z -= a.z*r;
      b = normalize(b);
    }
  }
  
  stl_position position(stl_vertex vertex, float tolerance = TOLERANCE) {
    float result = x*vertex.x + y*vertex.y + z*vertex.z + d;
    if (ABS(result) <= tolerance) return on;
    if (result > 0) return above;
    return below;
  }
  
  stl_vertex intersection(stl_vertex a, stl_vertex b) {
    stl_vertex ab; // vector from A to B
    ab.x = b.x-a.x;
    ab.y = b.y-a.y;
    ab.z = b.z-a.z;
    float t = - (a.x*x + a.y*y + a.z*z + d) / (ab.x*x + ab.y*y + ab.z*z);
    
    stl_vertex result;
    result.x = a.x + ab.x*t;
    result.y = a.y + ab.y*t;
    result.z = a.z + ab.z*t;
    return result;
  }
  
  stl_vertex to_2D(stl_vertex vertex, stl_vertex origin) {
    stl_vertex ov;
    ov.x = vertex.x-origin.x;
    ov.y = vertex.y-origin.y;
    ov.z = vertex.z-origin.z;
    
    stl_vertex result;
    result.x = scalar(a,ov);
    result.y = scalar(b,ov);
    result.z = 0;
    return result;
  }
  
  stl_vertex to_3D(stl_vertex vertex, stl_vertex origin) {
    stl_vertex result;
    
    result.x = origin.x + a.x*vertex.x + b.x*vertex.y;
    result.y = origin.y + a.y*vertex.x + b.y*vertex.y;
    result.z = origin.z + a.z*vertex.x + b.z*vertex.y;
    return result;
  }
};

struct stl_vertex_pair {
  stl_vertex x;
  stl_vertex y;
  stl_vertex_pair(stl_vertex x, stl_vertex y) {
    this->x=x;
    this->y=y;
  }
  bool operator<(const stl_vertex_pair& other) const  {
    if (x.x == other.x.x) {
      if (x.y == other.x.y) {
        if (x.z == other.x.z) {
          if (y.x == other.y.x) {
            if (y.y == other.y.y) {
              return y.z < other.y.z;
            } else return y.y < other.y.y;
          } else return y.x < other.y.x;
        } else return x.z < other.x.z;
      } else return x.y < other.x.y;
    } else return x.x < other.x.x;
  }
};

stl_facet semifacet(stl_facet original, stl_vertex a, stl_vertex b, stl_vertex c) {
  stl_facet f;
  f.vertex[0] = a;
  f.vertex[1] = b;
  f.vertex[2] = c;
  f.normal = original.normal;
  f.extra[0] = original.extra[0];
  f.extra[1] = original.extra[1];
  return f;
}

void simple_cut(stl_vertex zero, stl_vertex one, stl_vertex two, stl_facet facet, stl_plane plane,
              std::deque<stl_facet> &first, std::deque<stl_facet> &second,
              std::set<stl_vertex_pair> &border) {
  stl_vertex middle = plane.intersection(one, two);
  first.push_back(semifacet(facet, middle, zero, one));
  second.push_back(semifacet(facet, middle, two, zero));
  border.insert(stl_vertex_pair(one,middle));
}

void complex_cut(stl_vertex zero, stl_vertex one, stl_vertex two, stl_facet facet, stl_plane plane,
              std::deque<stl_facet> &first, std::deque<stl_facet> &second,
              std::set<stl_vertex_pair> &border) {
  stl_vertex one_middle = plane.intersection(zero, one);
  stl_vertex two_middle = plane.intersection(zero, two);
  first.push_back(semifacet(facet, zero, one_middle, two_middle));
  second.push_back(semifacet(facet, one_middle, one, two));
  second.push_back(semifacet(facet, one_middle, two, two_middle));
  border.insert(stl_vertex_pair(one_middle,two_middle));
}

void separate(stl_facet facet, stl_plane plane,
              std::deque<stl_facet> &upper, std::deque<stl_facet> &lower,
              std::set<stl_vertex_pair> &border) {
  stl_position pos[3];
  size_t aboves = 0;
  size_t belows = 0;
  size_t ons = 0;
  
  for (size_t i = 0; i < 3; i++) {
    pos[i] = plane.position(facet.vertex[i]);
    if (pos[i]==above) aboves++;
    else if (pos[i]==below) belows++;
    else ons++;
  }
    
  // All vertexes are above the plane
  if (aboves == 3) {
    upper.push_back(facet);
    return;
  }
  
  // All vertexes are below the plane
  if (belows == 3) {
    lower.push_back(facet);
    return;
  }
  
  // All vertexes are on the plane
  if (ons == 3) return;
  
  // 2 vertexes are on the plane
  if (ons == 2) {
    for (size_t i = 0; i < 3; i++) {
      if (pos[i] == above) {
        upper.push_back(facet);
        border.insert(stl_vertex_pair(facet.vertex[(i+1)%3],facet.vertex[(i+2)%3]));
      } else if (pos[i] == below) {
        lower.push_back(facet);
        border.insert(stl_vertex_pair(facet.vertex[(i+2)%3],facet.vertex[(i+1)%3]));
      }
    }
    return;
  }
  
  // 1 vertex is on the plane
  if (ons == 1) {
    stl_vertex zero, one, two;
    stl_position onepos;
    for (size_t i = 0; i < 3; i++) {
      if (pos[i] == on) {
        zero = facet.vertex[i];
        one = facet.vertex[(i+1)%3];
        onepos = pos[(i+1)%3];
        two = facet.vertex[(i+2)%3];
        break;
      }
    }
    if (aboves == 2) {
      upper.push_back(facet);
      return;
    }
    if (belows == 2) {
      lower.push_back(facet);
      return;
    }
    if (onepos == above)
      simple_cut(zero, one, two, facet, plane, upper, lower, border);
    else
      simple_cut(zero, one, two, facet, plane, lower, upper, border);
    return;
  }
  
  // no vertexes on the plane
  if  (aboves == 1) { // belows == 2
    for (size_t i = 0; i < 3; i++) {
      if (pos[i] == above) {
        complex_cut(facet.vertex[i], facet.vertex[(i+1)%3], facet.vertex[(i+2)%3], facet, plane, upper, lower, border);
        return;
      }
    }
  }
  // belows == 1, aboves == 2
  for (size_t i = 0; i < 3; i++) {
    if (pos[i] == below) {
      complex_cut(facet.vertex[i], facet.vertex[(i+1)%3], facet.vertex[(i+2)%3], facet, plane, lower, upper, border);
      return;
    }
  }
}

void export_stl(std::deque<stl_facet> facets, const char* name) {
  stl_file stl_out;
  stl_out.stats.type = inmemory;
  stl_out.stats.number_of_facets = facets.size();
  stl_out.stats.original_num_facets = stl_out.stats.number_of_facets;
  stl_out.v_indices = NULL;
  stl_out.v_shared = NULL;
  stl_out.neighbors_start = NULL;
  stl_clear_error(&stl_out);
  stl_allocate(&stl_out);
  
  int first = 1;
  for (std::deque<stl_facet>::const_iterator facet = facets.begin(); facet != facets.end(); facet++) {
    stl_out.facet_start[facet - facets.begin()] = *facet;
    stl_facet_stats(&stl_out, *facet, first);
    first = 0;
  }
  
  stl_write_ascii(&stl_out, name, "stlcut");
  stl_clear_error(&stl_out);
  stl_close(&stl_out);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " file.stl" << std::endl;
    return 1;
  }
  stl_file stl_in;
  stl_open(&stl_in, argv[1]);
  stl_exit_on_error(&stl_in);
  
  std::set<stl_vertex_pair> border;
  std::deque<stl_facet> upper, lower;
  
  for (int i = 0; i < stl_in.stats.number_of_facets; i++)
    separate(stl_in.facet_start[i], stl_plane(0,0,1,0), upper, lower, border);
  
  stl_close(&stl_in);
  
  export_stl(upper, "upper.stl");
  export_stl(lower, "lower.stl");
  
  return 0;
}

