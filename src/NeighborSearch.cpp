#include "..\headers\NeighborSearch.h"

// Standard Neighbor search
void QuadraticSearch() {
    for (int i = 0; i < particles.size(); i++) {
        particles[i].neighbors.clear();

        for (int j = 0; j < particles.size(); j++) {
            double distance = std::sqrt((particles[j].position - particles[i].position).squaredNorm());  

            if (distance < parameters.support) {
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
                        double distance = std::sqrt((p->position - particles[i].position).squaredNorm());

                        if (distance < parameters.support) {
							particles[i].neighbors.push_back(p);
						}
                    }
				}
            }
        }
	}
}

void NSghostsGrid() {
    // search neighbors for ghost particles only
    #pragma omp parallel for
    for (int i = 0; i < ghostParticles.size(); i++) {
        ghostParticles[i].neighbors.clear();
		const Eigen::Vector3i cellNumber = ghostParticles[i].getCellNumber();

        for (int j = cellNumber.x() - 1; j <= cellNumber.x() + 1; j++) {
            for (int k = cellNumber.y() - 1; k <= cellNumber.y() + 1; k++) {
                for (int l = cellNumber.z() - 1; l <= cellNumber.z() + 1; l++) {
					if (j < 0 || j >= GRID_WIDTH || 
                        k < 0 || k >= GRID_HEIGHT || 
						l < 0 || l >= GRID_DEPTH) { continue; }

                    for (auto& p : grid[j + k * GRID_WIDTH + l * GRID_WIDTH * GRID_HEIGHT].cellParticles) {
						double distance = std::sqrt((p->position - ghostParticles[i].position).squaredNorm());

                        if (distance < parameters.support) {
							ghostParticles[i].neighbors.push_back(p);
						}
					}
				}
			}
		}
    }
}

void NSUniformGrid2D() {
	GridUpdate2D();

	#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
		particles2D[i].neighbors.clear();
		const Eigen::Vector2i cellNumber = particles2D[i].getCellNumber();

        for (int j = cellNumber.x() - 1; j <= cellNumber.x() + 1; j++) {
            for (int k = cellNumber.y() - 1; k <= cellNumber.y() + 1; k++) {
				if (j < 0 || j >= GRID_WIDTH || k < 0 || k >= GRID_HEIGHT) { continue; }

                for (auto& p : grid2D[j + k * GRID_WIDTH].cellParticles) {
					double distance = std::sqrt((p->position - particles2D[i].position).squaredNorm());

                    if (distance < parameters.support) {
						particles2D[i].neighbors.push_back(p);
					}
				}
			}
		}
	}
}
