#include "..\headers\Grid.h"

Grid2D::Grid2D() :
	gridWidth(std::ceil(parameters.windowSize.width / parameters.spacing / 2)),
	gridHeight(std::ceil(parameters.windowSize.height / parameters.spacing / 2)),
	cellSize(parameters.spacing * 2) {
	initializeGrid();
}

Grid2D::Grid2D(double size) :
    gridWidth(std::ceil(parameters.windowSize.width / size)),
    gridHeight(std::ceil(parameters.windowSize.height / size)),
    cellSize(size) {
	initializeGrid();
}

void Grid2D::initializeGrid() {
	std::cout << "Using uniform grid with " << gridWidth << "x" << gridHeight << " cells" << std::endl;
    cells2D.resize(gridWidth * gridHeight);
}

void Grid2D::updateGrid() {
	for (auto& cell : cells2D) {
		cell.cellParticles.clear();
	}
    for (auto& p : particles2D) {
		const Eigen::Vector2i cellNumber = (p.position / cellSize).cast<int>();
        if (cellNumber.x() < 0 || cellNumber.x() >= gridWidth ||
            cellNumber.y() < 0 || cellNumber.y() >= gridHeight) {
            continue;
        }
        cells2D[cellNumber.x() + cellNumber.y() * gridWidth].cellParticles.push_back(&p);
    }
}

void Grid2D::neighborSearch(std::vector<Particle2D>& particles2D) {
#pragma omp parallel for
    for (int i = 0; i < particles2D.size(); i++) {
        particles2D[i].neighbors.clear();
        const Eigen::Vector2i cellNumber = (particles2D[i].position / cellSize).cast<int>();

        for (int j = cellNumber.x() - 1; j <= cellNumber.x() + 1; j++) {
            for (int k = cellNumber.y() - 1; k <= cellNumber.y() + 1; k++) {
                if (j < 0 || j >= gridWidth || k < 0 || k >= gridHeight) { continue; }

                for (auto& p : cells2D[j + k * gridWidth].cellParticles) {
                    double distance = std::sqrt((p->position - particles2D[i].position).squaredNorm());

                    if (distance < parameters.support) {
                        particles2D[i].neighbors.push_back(p);
                    }
                }
            }
        }
    }
}

Grid::Grid(double size)
    : Grid2D(size), gridDepth(std::ceil(parameters.windowSize.depth / size)) {
    initializeGrid();
}

void Grid::initializeGrid() {
    std::cout << "Using uniform grid with " << gridWidth << "x" << gridHeight << "x" << gridDepth << " cells" 
        << std::endl;
    cells.resize(gridWidth * gridHeight * gridDepth);
}

void Grid::updateGrid() {
    for (auto& cell : cells) {
		cell.cellParticles.clear();
    }
    for (auto& p : particles) {
        const Eigen::Vector3i cellNumber = (p.position / cellSize).cast<int>();
        if (cellNumber.x() < 0 || cellNumber.x() >= gridWidth ||
            cellNumber.y() < 0 || cellNumber.y() >= gridHeight ||
            cellNumber.z() < 0 || cellNumber.z() >= gridDepth) {
			continue;
		}
		cells[cellNumber.x() + cellNumber.y() * gridWidth +
            cellNumber.z() * gridWidth * gridHeight].cellParticles.push_back(&p);
	}
}

void Grid::neighborSearch(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (int i = 0; i < particles.size(); i++) {
        particles[i].neighbors.clear();
        const Eigen::Vector3i cellNumber = (particles[i].position / cellSize).cast<int>();

        for (int j = cellNumber.x() - 1; j <= cellNumber.x() + 1; j++) {
            for (int k = cellNumber.y() - 1; k <= cellNumber.y() + 1; k++) {
                for (int l = cellNumber.z() - 1; l <= cellNumber.z() + 1; l++) {
                    if (j < 0 || j >= gridWidth || k < 0 || k >= gridHeight || l < 0 || l >= gridDepth) { continue; }

                    for (auto& p : cells[j + k * gridWidth + l * gridWidth * gridHeight].cellParticles) {
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

// Initialize uniformed grid of fluid particles
void InitFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            for (int k = 0; k < parameters.particlesPerDimension.z; k++) {
                Particle p;

                p.position = Eigen::Vector3d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing, (k + 2) * parameters.spacing);
				p.ID = particles.size();
                particles.push_back(p);
            }
        }
    }
}

void InitGhostFluid() {
    for (int i = 0; i < parameters.particlesPerDimension.x; i++) {
        for (int j = 0; j < parameters.particlesPerDimension.y; j++) {
            Particle p;

            p.position = Eigen::Vector3d((i + 2) * parameters.spacing, (j + 2) * parameters.spacing, parameters.windowSize.depth / 2);
            ghostParticles.push_back(p);
        }
    }
}

// Initialize boundaries
void InitBoundaries() {
    int width = parameters.windowSize.width / parameters.spacing - 1;
    int hight = parameters.windowSize.height / parameters.spacing - 1;
    int depth = (parameters.windowSize.depth / parameters.spacing - 1) / 2;

    for (float i = 0; i < width; i += 0.5) {
        for (float j = 0; j < hight; j += 0.5) {
            for (float k = 0; k < depth; k += 0.5) {
                if (i < 0.5 || i > width - 1 || j < 0.5 || j > hight - 1 || k < 0.5 || k > depth - 1) {
                    Particle p;

                    p.position = Eigen::Vector3d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing, (k + 1) * parameters.spacing);
                    p.isFluid = false;
                    p.ID = particles.size();
                    parameters.boundaryTestID = p.ID;

                    particles.push_back(p);
                }
			}
        }
    }

    /*for (float i = 0; i < width; i ++) {
        for (float j = 0; j < hight; j ++) {
            for (float k = 0; k < depth; k ++) {
                if (i < 1 || i > width - 2 || j < 1 || j > hight - 2 || k < 1 || k > depth - 2) {
                    Particle p;

                    p.position = Eigen::Vector3d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing, (k + 1) * parameters.spacing);
                    p.isFluid = false;
                    p.ID = particles.size();
                    particles.push_back(p);
                }
            }
        }
    }*/
}

void MovingBoundary() {
    for (int i = -20; i < 0; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = 0; k < (parameters.windowSize.depth / parameters.spacing - 1) / 2; k++) {
				Particle p;

				p.position = Eigen::Vector3d(i * parameters.spacing, j * parameters.spacing, (k + 1) * parameters.spacing);
				p.isFluid = false;
				p.ID = particles.size();
				p.velocity = Eigen::Vector3d(20, 0.0, 0.0);

				particles.push_back(p);
			}
        }
    }
}

void InitFluidForRotatingTest2D() {
    for (int i = 1; i < parameters.windowSize.width / parameters.spacing / 4; i++) {
        for (int j = (parameters.windowSize.height / parameters.spacing - 1) / 4 + 1; j < (parameters.windowSize.height / parameters.spacing - 1) / 4 + 35; j++) {
            Particle2D p;

            p.position = Eigen::Vector2d((i + 1) * parameters.spacing, (j + 1) * parameters.spacing);
            p.ID = particles2D.size();

            particles2D.push_back(p);
        }
    }
}
