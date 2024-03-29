#include "..\headers\Grid.h"

std::vector<Particle> particles;
std::vector<GridCell> grid;

int GRID_WIDTH;
int GRID_HEIGHT;
int GRID_DEPTH;

// Initialize uniformed grid of fluid particles
void InitFluid(const int l) {
    int depth = SCENE_DEPTH / SPACING - 1;

    for (int i = 0; i < PARTICLES_PER_DIMENSION * l; i++) {
        for (int j = 0; j < PARTICLES_PER_DIMENSION; j++) {
            for (int k = 3; k <= depth - 4; k++) {
                Particle p;

				p.position = Eigen::Vector3f((i + 4) * SPACING, (j + 4) * SPACING, (k + 1) * SPACING);
				p.ID = particles.size();
				particles.push_back(p);
            }
        }
    }
}

// Initialize boundaries
void InitBoundaries() {
    int width = WINDOW_WIDTH / SPACING - 1;
    int hight = WINDOW_HEIGHT / SPACING - 1;
    int depth = SCENE_DEPTH / SPACING - 1;

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < hight; j++) {
            for (int k = 0; k < depth; k++) {
                if (i < 3 || i > width - 4 || j < 3 || j > hight - 4 || k < 3 || k > depth - 4) {
                    Particle p;

                    p.position = Eigen::Vector3f((i + 1) * SPACING, (j + 1) * SPACING, (k + 1) * SPACING);
                    p.isFluid = false;
                    p.ID = particles.size();

                    particles.push_back(p);
                }
			}
        }
    }
}

void UniformGrid() {
    GRID_WIDTH = std::ceil(WINDOW_WIDTH / CELL_SIZE);
    GRID_HEIGHT = std::ceil(WINDOW_HEIGHT / CELL_SIZE);
    GRID_DEPTH = std::ceil(SCENE_DEPTH / CELL_SIZE);

    std::cout << "Using uniform grid with " << GRID_WIDTH  << "x" 
                                            << GRID_HEIGHT << "x" 
                                            << GRID_DEPTH  << " cells" << std::endl;

    grid.resize(GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH);
}

void GridUpdate() {
    // Clear grid
    for (int i = 0; i < grid.size(); i++) {
        grid[i].cellParticles.clear();
    }

    for (auto& p : particles) {
        const Eigen::Vector3i cellNumber = p.getCellNumber();
        if (cellNumber.x() < 0 || cellNumber.x() >= GRID_WIDTH || 
            cellNumber.y() < 0 || cellNumber.y() >= GRID_HEIGHT || 
            cellNumber.z() < 0 || cellNumber.z() >= GRID_DEPTH) {
            continue;
        }
        grid[cellNumber.x() + cellNumber.y() * GRID_WIDTH + 
             cellNumber.z() * GRID_WIDTH * GRID_HEIGHT].cellParticles.push_back(&p);
    }
}
