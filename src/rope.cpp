#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass,
           float k, vector<int> pinned_nodes) {
  // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and
  // containing `num_nodes` nodes.
  //
  // New masses can be created (see mass.h):
  // Mass *m = new Mass(position, mass, bool)
  // and should be added to the masses vector.
  //
  // Springs can be created (see spring.h):
  // Spring *s = new Spring(mass_a, mass_b, spring_const)
  // and should be added to the springs vector.
  //
  // Masses corresponding to indices in pinned_nodes
  // should have their pinned field set to true.
    
    // You will have two vectors to store mass and spring.
    // i.e. mass.push_back(mass);
    //i.e. springs.push_back(spring)
    //Start by computing the offset between nodes
    Vector2D offset = (end-start) / (num_nodes - 1);
    
    for (int i = 0; i < num_nodes; i++) {
        // do something for mass
        Vector2D position = start + (i * offset);
        Mass *m = new Mass(position, node_mass, false);
        masses.push_back(m);
        
        if(i > 0) {
            // do something for spring
            Spring *s = new Spring(m, masses[i-1], k);
            springs.push_back(s);
        }
    }
    //final cehck: check if i is a pinned node
    for (auto &i : pinned_nodes) {
        masses[i]->pinned = true;
    }

}

void Rope::simulateEuler(float delta_t, Vector2D gravity) {
  for (auto &s : springs) {
    // TODO (Part 2.1): Use Hooke's law to calculate the force on a node
      // s->m2-> position, s->m1->position
      Vector2D b = s->m2-> position;
      Vector2D a = s->m1->position;
      
      Vector2D diff = b - a;
      
      Vector2D spring_force = s->k * diff.unit() * (diff.norm() - s->rest_length);
      s->m1->forces += spring_force;
      s->m2->forces -= spring_force;
      
    // TODO (Part 4.1): Add damping forces
      Vector2D relative_velocity= s->m2->velocity - s->m1->velocity;
      Vector2D damping_force = 0.5 * relative_velocity;
      s->m1->forces += damping_force;
      s->m2->forces -= damping_force;

    // TODO Apply forces as appropriate.

  }

  for (auto &m : masses) {
    if (!m->pinned) {
      // TODO (Part 2.1): Add the force due to gravity, then compute the new
      // velocity and position
        m->velocity += ((m->forces + gravity)/m->mass) * delta_t;
        m->position += m->velocity * delta_t;
    }

    // TODO Reset all forces on each mass
      m->forces = Vector2D(0,0);

  }
}

void Rope::simulateVerlet(float delta_t, Vector2D gravity) {
  // TODO (Part 3.1): Clear forces
  for (auto &s : springs) {
    // TODO (Part 3.1): Simulate one timestep of the rope using explicit Verlet
    
  }

  for (auto &m : masses) {
    if (!m->pinned) {
      Vector2D temp_position = m->position;

      // TODO (Part 3.1): Set the new position of the rope mass
      // TODO (Part 4.2): Add global Verlet damping

    }
  }
}
}
