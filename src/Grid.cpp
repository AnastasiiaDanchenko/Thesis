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

Grid::Grid() {
	gridDepth = std::ceil(parameters.windowSize.depth / parameters.spacing / 2);
	initializeGrid();
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
