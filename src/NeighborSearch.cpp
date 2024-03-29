#include "..\headers\NeighborSearch.h"

// Standard Neighbor search
void QuadraticSearch() {
    for (int i = 0; i < particles.size(); i++) {
        particles[i].neighbors.clear();

        for (int j = 0; j < particles.size(); j++) {
            float distance = std::sqrt((particles[j].position - particles[i].position).squaredNorm());  

            if (distance < SUPPORT) {
                particles[i].neighbors.push_back(&particles[j]);
            }
        }
    }
}

void NSUniformGrid() {
    GridUpdate();

	#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
		particles[i].neighbors.clear();
        const Eigen::Vector3i cellNumber = particles[i].getCellNumber();

        for (int j = cellNumber.x() - 1; j <= cellNumber.x() + 1; j++) {
            for (int k = cellNumber.y() - 1; k <= cellNumber.y() + 1; k++) {
                for (int l = cellNumber.z() - 1; l <= cellNumber.z() + 1; l++) {
                    if (j < 0 || j >= GRID_WIDTH || 
                        k < 0 || k >= GRID_HEIGHT || 
                        l < 0 || l >= GRID_DEPTH) { continue; }

                    for (auto& p : grid[j + k * GRID_WIDTH + l * GRID_WIDTH * GRID_HEIGHT].cellParticles) {
                        float distance = std::sqrt((p->position - particles[i].position).squaredNorm());

                        if (distance < SUPPORT) {
							particles[i].neighbors.push_back(p);
						}
                    }
				}
            }
        }
	}
}
