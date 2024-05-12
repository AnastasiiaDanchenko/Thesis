#include "..\headers\Grid.h"

std::vector<Particle> particles;
std::vector<Particle2D> particles2D;

std::vector<GridCell> grid;
std::vector<GridCell2D> grid2D;

int GRID_WIDTH;
int GRID_HEIGHT;
int GRID_DEPTH;

// Initialize uniformed grid of fluid particles
void InitFluid() {
    int depth = SCENE_DEPTH / SPACING - 1;

    for (int i = 0; i < PARTICLES_X; i++) {
        for (int j = 0; j < PARTICLES_Y; j++) {
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

void InitFluid2D() {
	int depth = SCENE_DEPTH / SPACING - 1;

    for (int i = 0; i < PARTICLES_X; i++) {
        for (int j = 0; j < PARTICLES_Y; j++) {
			Particle2D p;

			p.position = Eigen::Vector2f((i + 4) * SPACING, (j + 4) * SPACING);
			p.ID = particles2D.size();
            particles2D.push_back(p);
		}
	}
}

void InitBoundaries2D() {
	int width = (WINDOW_WIDTH / 2) / SPACING - 1;
	int hight = (WINDOW_HEIGHT / 2) / SPACING - 1;

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < hight; j++) {
            if (i < 3 || i > width - 4 || j < 3 || j > hight - 4) {
				Particle2D p;

				p.position = Eigen::Vector2f((i + 1) * SPACING, (j + 1) * SPACING);
				p.isFluid = false;
				p.ID = particles2D.size();

                particles2D.push_back(p);
			}
		}
	}
}

void UniformGrid2D() {
	GRID_WIDTH = std::ceil(WINDOW_WIDTH / CELL_SIZE);
	GRID_HEIGHT = std::ceil(WINDOW_HEIGHT / CELL_SIZE);

	std::cout << "Using uniform grid with " << GRID_WIDTH << "x"
		<< GRID_HEIGHT << " cells" << std::endl;

	grid2D.resize(GRID_WIDTH * GRID_HEIGHT);
}

void GridUpdate2D() {
	// Clear grid
    for (int i = 0; i < grid2D.size(); i++) {
		grid2D[i].cellParticles.clear();
	}

    for (auto& p : particles2D) {
		const Eigen::Vector2i cellNumber = p.getCellNumber();
        if (cellNumber.x() < 0 || cellNumber.x() >= GRID_WIDTH ||
            cellNumber.y() < 0 || cellNumber.y() >= GRID_HEIGHT) {
			continue;
		}
		grid2D[cellNumber.x() + cellNumber.y() * GRID_WIDTH].cellParticles.push_back(&p);
	}
}
