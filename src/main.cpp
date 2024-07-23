#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <boost/polygon/voronoi.hpp>
#include <SFML/Graphics.hpp>
#include <unordered_map>
#include <stack>
#include <cmath>
#include <unordered_set>
#include <algorithm>
#include <numeric>


using namespace boost::polygon;
typedef double coordinate_type;



struct Point {
    coordinate_type x, y;
    Point(coordinate_type x_, coordinate_type y_) : x(x_), y(y_) {}
    
};

struct Cluster {
    std::vector<Point> points;
    std::vector<int> indices;
    Point centroid;
    double radius;

    Cluster() : centroid(Point(0, 0)), radius(0) {}
};

using ClusterDictionary = std::unordered_map<double, Cluster>;

namespace boost {
namespace polygon {

template <>
struct geometry_concept<Point> {
    typedef point_concept type;
};

template <>
struct point_traits<Point> {
    typedef double coordinate_type;

    static inline coordinate_type get(const Point& point, orientation_2d orient) {
        return (orient == HORIZONTAL) ? point.x : point.y;
    }
};

}  // namespace polygon
}  // namespace boost

void buildVoronoiDiagram(const std::vector<Point>& points, voronoi_diagram<double>& vd) {
    construct_voronoi(points.begin(), points.end(), &vd);
}

// Función para calcular la distancia euclidiana entre dos puntos
double distance(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
}

double poisson_cdf(double k, double lambda) {
    double cdf = 0.0;
    for (int i = 0; i <= k; ++i) {
        cdf += (std::pow(lambda, i) * std::exp(-lambda)) / tgamma(i + 1);
    }
    return cdf;
}

std::vector<Point> readPointsFromCSV(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening CSV file: " << filename << std::endl;
        return points;
    }

    std::string line;
    bool is_header = true;
    while (std::getline(file, line)) {
        if (is_header) {
            is_header = false;
            continue; // Skip header line
        }

        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(iss, token, ',')) {
            tokens.push_back(token);
        }
        if (tokens.size() != 2) {
            std::cerr << "Error: invalid line in CSV file: " << line << std::endl;
            continue;
        }

        try {
            double x = std::stod(tokens[0]);
            double y = std::stod(tokens[1]);
            points.emplace_back(x, y);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error converting values in CSV file: " << e.what() << std::endl;
        }
    }
    return points;
}

void color_infection(voronoi_diagram<double>::cell_type* cell, const voronoi_diagram<double>& vd, double color, const std::vector<Point>& points, double max_distance) {
    cell->color(color);

    const voronoi_edge<double>* edge = cell->incident_edge();
    const voronoi_edge<double>* start = edge;
    do {
        const voronoi_diagram<double>::cell_type* neighbor = edge->twin()->cell();
        if (neighbor->color() == 1) {
            std::size_t neighbor_index = neighbor->source_index();
            double dist = distance(points[cell->source_index()], points[neighbor_index]);
            if (dist <= max_distance) {
                color_infection(const_cast<voronoi_diagram<double>::cell_type*>(neighbor), vd, color, points, max_distance);
            }
        }
        edge = edge->next();
    } while (edge != start);
}

double calculateVoronoiCellArea(const voronoi_diagram<double>::cell_type& cell, const voronoi_diagram<double>& vd) {
    if (!cell.contains_point() && !cell.contains_segment()) {
        return 0.0; // Invalid cell
    }
    double area = 0.0;
    const voronoi_edge<double>* edge = cell.incident_edge();
    const voronoi_edge<double>* start = edge;
    do {
        if (!edge) {
            return 0.0; // Handle null pointer
        }
        
        if (edge->is_primary()) {
            const voronoi_vertex<double>* v0 = edge->vertex0();
            const voronoi_vertex<double>* v1 = edge->vertex1();
            if (v0 && v1) {
                area += (v0->x() * v1->y() - v1->x() * v0->y());
            }
        }
        edge = edge->next();
    } while (edge != start);
    return std::abs(area) / 2.0;
}

ClusterDictionary collectClusters(const voronoi_diagram<double>& vd, const std::vector<Point>& points) {
    ClusterDictionary clusters;

    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        int color = cell_it->color();
        if (color > 1) { // Ignorar celdas no coloreadas o aisladas
            std::size_t index = cell_it->source_index();
            clusters[color].points.push_back(points[index]);
            clusters[color].indices.push_back(cell_it - vd.cells().begin()); // Guardar el índice de la celda
        }
    }

    return clusters;
}



Point calculateCentroid(const std::vector<Point>& points) {
    double sum_x = 0;
    double sum_y = 0;

    for (const auto& point : points) {
        sum_x += point.x;
        sum_y += point.y;
    }

    return Point(sum_x / points.size(), sum_y / points.size());
}

double calculateRadius(const Point& centroid, const std::vector<Point>& points) {
    double max_distance = 0;

    for (const auto& point : points) {
        double distance = std::sqrt(std::pow(centroid.x - point.x, 2) + std::pow(centroid.y - point.y, 2));
        if (distance > max_distance) {
            max_distance = distance;
        }
    }

    return max_distance;
}

void calculateClusterProperties(voronoi_diagram<double>& vd, const std::vector<Point>& points, const std::string& output_file) {
    auto clusters = collectClusters(vd, points);

    for (auto& [color, cluster] : clusters) {
        cluster.centroid = calculateCentroid(cluster.points);
        cluster.radius = calculateRadius(cluster.centroid, cluster.points);
    }

    // Guardar los resultados en un archivo CSV
    std::ofstream file(output_file);
    file << "ClusterID,CentroidX,CentroidY,Radius\n";
    for (const auto& [color, cluster] : clusters) {
        file << color << "," << cluster.centroid.x << "," << cluster.centroid.y << "," << cluster.radius << "\n";
    }
    file.close();
}

void colorVoronoiCells(voronoi_diagram<double>& vd, double area_threshold, double xmin, double xmax, double ymin, double ymax,
 double model_distance, const std::vector<Point>& points) {
    std::set<std::size_t> small_cells;

    // Identificar celdas pequeñas
    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        if (cell_it->contains_point()) {
            bool is_external = false;
            const voronoi_edge<double>* edge = cell_it->incident_edge();
            const voronoi_edge<double>* start = edge;
            do {
                const voronoi_vertex<double>* v = edge->vertex0();
                if (v && (v->x() < xmin || v->x() > xmax || v->y() < ymin || v->y() > ymax)) {
                    is_external = true;
                    break;
                }
                edge = edge->next();
            } while (edge != start);

            if (!is_external) {
                double area = calculateVoronoiCellArea(*cell_it, vd);
                std::size_t index = cell_it->source_index();
                if (area < area_threshold) {
                    small_cells.insert(cell_it - vd.cells().begin());
                    cell_it->color(1); // Marcar celdas pequeñas con color 1 (rojo)
                }
            }
        }
    }

    // Eliminar celdas pequeñas aisladas
    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        if (cell_it->color() == 1) {
            int neighboring_colored_cells = 0;
            const voronoi_edge<double>* edge = cell_it->incident_edge();
            const voronoi_edge<double>* start = edge;
            do {
                if (edge->is_primary()) {
                    const voronoi_vertex<double>* v = edge->vertex0();
                    if (v) {
                        const voronoi_cell<double>* neighbor_cell = edge->twin()->cell();
                        if (neighbor_cell && neighbor_cell->color() == 1) {
                            ++neighboring_colored_cells;
                        }
                    }
                }
                edge = edge->next();
            } while (edge != start);

            if (neighboring_colored_cells < 2) {
                cell_it->color(0); // Despintar celdas aisladas
                small_cells.erase(cell_it - vd.cells().begin());
            }
        }
    }

    // Identificar celdas vecinas a las del cluster
    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        if (cell_it->color() != 1) {
            double min_distance = std::numeric_limits<double>::max();
            std::size_t cell_index = cell_it->source_index();
            const Point& cell_point = points[cell_index];

            // Recorrer celdas vecinas
            const voronoi_edge<double>* edge = cell_it->incident_edge();
            const voronoi_edge<double>* start = edge;
            do {
                if (edge->is_primary()) {
                    const voronoi_vertex<double>* v0 = edge->vertex0();
                    if (v0) {
                        std::size_t neighbor_index = v0->incident_edge()->cell()->source_index();
                        if (small_cells.find(neighbor_index) != small_cells.end()) {
                            const Point& neighbor_point = points[neighbor_index];
                            double dist = distance(cell_point, neighbor_point);
                            if (dist < min_distance) {
                                min_distance = dist;
                            }
                        }
                    }
                }
                edge = edge->next();
            } while (edge != start);

            if (min_distance <= model_distance) {
                cell_it->color(2); // Marcar celdas vecinas con color 2 (azul)
            }
        }
    }
}

sf::Color getColorForCluster(std::size_t index) {
    // Define the maximum number of colors and the base hue range
    const std::size_t maxColors = 100;
    float hue = static_cast<float>(index % maxColors) / static_cast<float>(maxColors);
    
    // Convert hue to RGB color
    sf::Uint8 r, g, b;
    float h = hue * 6.0f;
    float f = h - static_cast<int>(h);
    float p = 0.0f;
    float q = 1.0f - f;
    float t = f;
    
    switch (static_cast<int>(h)) {
        case 0: r = 255; g = static_cast<sf::Uint8>(t * 255); b = static_cast<sf::Uint8>(p * 255); break;
        case 1: r = static_cast<sf::Uint8>(q * 255); g = 255; b = static_cast<sf::Uint8>(p * 255); break;
        case 2: r = static_cast<sf::Uint8>(p * 255); g = 255; b = static_cast<sf::Uint8>(t * 255); break;
        case 3: r = static_cast<sf::Uint8>(p * 255); g = static_cast<sf::Uint8>(q * 255); b = 255; break;
        case 4: r = static_cast<sf::Uint8>(t * 255); g = static_cast<sf::Uint8>(p * 255); b = 255; break;
        case 5: r = 255; g = static_cast<sf::Uint8>(p * 255); b = static_cast<sf::Uint8>(q * 255); break;
    }
    
    return sf::Color(r, g, b, 128); // Use alpha for transparency
}

void visualizeVoronoiDiagram(const voronoi_diagram<double>& vd, const std::vector<Point>& points) {
    sf::RenderWindow window(sf::VideoMode(800, 600), "Voronoi Diagram");
    sf::View view = window.getView();
    const float moveSpeed = 20.0f; // Speed at which the view moves
    bool dragging = false;
    sf::Vector2i oldMousePos;
    while (window.isOpen()) {
        sf::Event event;
        
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Add) { // '+' key for zoom in
                    view.zoom(0.9f); // Zoom in by reducing the view size
                } else if (event.key.code == sf::Keyboard::Subtract) { // '-' key for zoom out
                    view.zoom(1.1f); // Zoom out by increasing the view size
                } else if (event.key.code == sf::Keyboard::Up) { // Up arrow key to move view up
                    view.move(0, -moveSpeed);
                } else if (event.key.code == sf::Keyboard::Down) { // Down arrow key to move view down
                    view.move(0, moveSpeed);
                } else if (event.key.code == sf::Keyboard::Left) { // Left arrow key to move view left
                    view.move(-moveSpeed, 0);
                } else if (event.key.code == sf::Keyboard::Right) { // Right arrow key to move view right
                    view.move(moveSpeed, 0);
                }
            }

            if (event.type == sf::Event::MouseWheelMoved) {
                if (event.mouseWheel.delta > 0) {
                    view.zoom(0.9f);
                } else {
                    view.zoom(1.1f);
                }
            }

            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    // Iniciar el arrastre
                    dragging = true;
                    oldMousePos = sf::Mouse::getPosition(window);
                }
            }

            if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    // Terminar el arrastre
                    dragging = false;
                }
            }

            if (event.type == sf::Event::MouseMoved) {
                if (dragging) {
                    sf::Vector2i newMousePos = sf::Mouse::getPosition(window);
                    sf::Vector2f oldWorldPos = window.mapPixelToCoords(oldMousePos);
                    sf::Vector2f newWorldPos = window.mapPixelToCoords(newMousePos);

                    sf::Vector2f delta = oldWorldPos - newWorldPos;
                    view.move(delta);

                    oldMousePos = newMousePos; // Actualizar la posición del mouse
                }
            }
        }

        window.clear(sf::Color::White);
        window.setView(view);

        // Draw the points
        for (const auto& point : points) {
        sf::ConvexShape triangle;
        int size = 4;  // Tamaño del triángulo
        triangle.setPointCount(3);  // Un triángulo tiene 3 puntos
        triangle.setPoint(0, sf::Vector2f(point.x, point.y - size));  // Punto superior
        triangle.setPoint(1, sf::Vector2f(point.x - size, point.y + size));  // Punto inferior izquierdo
        triangle.setPoint(2, sf::Vector2f(point.x + size, point.y + size));  // Punto inferior derecho
        triangle.setFillColor(sf::Color::Red);
        window.draw(triangle);
    }


        // Draw the Voronoi edges and highlight cluster cells
        for (auto edge_it = vd.edges().begin(); edge_it != vd.edges().end(); ++edge_it) {
            if (edge_it->is_infinite()) {
                continue;
            }
            const voronoi_vertex<double>* v0 = edge_it->vertex0();
            const voronoi_vertex<double>* v1 = edge_it->vertex1();
            if (v0 != nullptr && v1 != nullptr) {
                sf::Vertex line[] =
                {
                    sf::Vertex(sf::Vector2f(v0->x(), v0->y()), sf::Color::Black),
                    sf::Vertex(sf::Vector2f(v1->x(), v1->y()), sf::Color::Black)
                };
                window.draw(line, 2, sf::Lines);
            }
        }

        // Highlight colored cells
        for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
            std::size_t index = cell_it->source_index();
            const voronoi_edge<double>* edge = cell_it->incident_edge();
            const voronoi_edge<double>* start = edge;
            sf::Color cell_color;
            sf::ConvexShape convex;
            if (cell_it->color() == 1) {
                cell_color = sf::Color::Red;
            } else if (cell_it->color() > 1) {
                cell_color = sf::Color::Blue;
            } else {
                continue;
            }

            sf::ConvexShape polygon;
            if (cell_it->color() > 0) {
                cell_color = getColorForCluster(cell_it->color() - 1); // Adjust color index based on internal color
            }

            polygon.setFillColor(cell_color);
            
            std::vector<sf::Vector2f> vertices;

            do {
                const voronoi_vertex<double>* v = edge->vertex0();
                if (v != nullptr) {
                    vertices.push_back(sf::Vector2f(v->x(), v->y()));
                }
                edge = edge->next();
            } while (edge != start);

            polygon.setPointCount(vertices.size());
            for (size_t i = 0; i < vertices.size(); ++i) {
                polygon.setPoint(i, vertices[i]);
            }
            window.draw(polygon);
        }

        window.display();
    }
}







// Función para identificar celdas pequeñas
std::vector<std::size_t> findSmallCells(const voronoi_diagram<double>& vd, double threshold, double xmin, double xmax, double ymin, double ymax) {
    std::vector<std::size_t> small_cells;
    std::vector<double> areas;



    // Identificar y marcar las celdas densas
    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        if (cell_it->contains_point()) {
            bool is_external = false;
            const voronoi_edge<double>* edge = cell_it->incident_edge();
            const voronoi_edge<double>* start = edge;
            do {
                const voronoi_vertex<double>* v = edge->vertex0();
                if (v && (v->x() < xmin || v->x() > xmax || v->y() < ymin || v->y() > ymax)) {
                    is_external = true;
                    break;
                }
                edge = edge->next();
            } while (edge != start);

            if (!is_external) {
                double area = calculateVoronoiCellArea(*cell_it, vd);
                std::size_t index = cell_it->source_index();
                if (area < threshold) {
                    small_cells.push_back(cell_it - vd.cells().begin());
                    cell_it->color(1); // Marcar celdas densas con color 1 (rojo)
                } else {
                    cell_it->color(0); // Marcar celdas no densas con color 0
                }
            } else {
                cell_it->color(0); // Marcar celdas externas con color 0
            }
        } else {
            cell_it->color(0); // Marcar celdas no válidas con color 0
        }
    }

    return small_cells;
}

// Función para eliminar celdas pequeñas aisladas
void removeIsolatedSmallCells(voronoi_diagram<double>& vd, ClusterDictionary& clusters, int min_points_threshold) {
    std::vector<int> clusters_to_remove;

    // Verificar cada cluster en el mapa
    for (const auto& [color, cluster] : clusters) {
        if (cluster.points.size() < min_points_threshold) {
            clusters_to_remove.push_back(color);
        }
    }

    // Eliminar clusters marcados
    for (const double color : clusters_to_remove) {
        const Cluster& cluster = clusters[color];
        for (const int& index : cluster.indices) {
            // Encontrar la celda correspondiente en el diagrama de Voronoi
            vd.cells()[index].color(0); // Despintar la celda
        }
        clusters.erase(color); // Eliminar el cluster del mapa
    }
}



// Función para pintar celdas con colores distintos (etiquetar)
void labelClusters(voronoi_diagram<double>& vd, const std::vector<Point>& points, double max_distance) {
    double color_id = 2; // Comenzamos desde 2, ya que 1 es para celdas pequeñas y 0 es para celdas no-cluster
    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        // Nos aseguramos de que solo procesamos celdas con color 1
        if (cell_it->color() == 1) {
            // Propagamos el color a través de las celdas vecinas
            color_infection(const_cast<voronoi_diagram<double>::cell_type*>(&(*cell_it)), vd, color_id++, points, max_distance);
        }
    }
}

// Función para conectar celdas cercanas a clusters
void connectNearbyCells(voronoi_diagram<double>& vd, const std::vector<Point>& points, double model_distance, double xmin, double xmax, double ymin, double ymax) {

    

    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        std::size_t index = cell_it->source_index();
        if (cell_it->color() == 0 && points[index].x > xmin && points[index].x < xmax && points[index].y > ymin && points[index].y < ymax) {
            
            Point cell_center = points[cell_it->source_index()];
            double min_distance = std::numeric_limits<double>::max();
            double closest_color = 0;

            const voronoi_edge<double>* edge = cell_it->incident_edge();
            const voronoi_edge<double>* start = edge;

            do {
                const voronoi_vertex<double>* v0 = edge->vertex0();
                if (v0) {
                    const voronoi_diagram<double>::cell_type* neighbor = edge->twin()->cell();
                    std::size_t neighbor_index = neighbor->source_index();
                    
                    if (neighbor && neighbor->color() != 0) {
                        Point neighbor_center = points[neighbor_index];
                        double dist = distance(cell_center, neighbor_center);
                        if (dist < min_distance) {
                            min_distance = dist;
                            closest_color = neighbor->color();
                        }
                    }
                }
                
                edge = edge->next();
            } while (edge != start);
            
            if (min_distance <= model_distance * 0.6) {
                cell_it->color(closest_color);
            }
        }
    }
}

// Función principal para el proceso de clusterización
//void clusterizeVoronoiDiagram(voronoi_diagram<double>& vd, const std::vector<Point>& points, double xmin, double xmax, double ymin, double ymax,
// double model_distance, double size_threshold, ) {


//    removeIsolatedSmallCells(vd, );
//   labelClusters(vd, points, model_distance);
//    connectNearbyCells(vd, points, model_distance, xmin,  xmax,  ymin, ymax);
//}

double dp(double a, double coefficient) {
        return coefficient * std::pow(a, 3) * std::exp(-4 * a);
}

void saveToCSV(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename) {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << "Area/Promedio,dp(a)\n";
            for (size_t i = 0; i < x.size(); ++i) {
                file << x[i] << "," << y[i] << "\n";
            }
            file.close();
        }
    }

int main() {
    std::string filename = "points.csv"; // Replace with your CSV filename
    std::vector<Point> points = readPointsFromCSV(filename);

    voronoi_diagram<double> vd;
    buildVoronoiDiagram(points, vd);

    // Calcula los límites de la mayoría de los puntos
    double xmin = points[0].x, xmax = points[0].x;
    double ymin = points[0].y, ymax = points[0].y;
    for (const auto& point : points) {
        xmin = std::min(xmin, point.x);
        xmax = std::max(xmax, point.x);
        ymin = std::min(ymin, point.y);
        ymax = std::max(ymax, point.y);
    }

    // Margen adicional para incluir todos los puntos internos
    double margin = -(xmax-xmin + ymax-ymin)/2 * 0.02;
    std::cout << "Margin: " << margin << std::endl;
    xmin -= margin;
    xmax += margin;
    ymin -= margin;
    ymax += margin;

    std::vector<double> areas;
    for (auto cell_it = vd.cells().begin(); cell_it != vd.cells().end(); ++cell_it) {
        if (cell_it->contains_point()) {
            bool is_external = false;
            const voronoi_edge<double>* edge = cell_it->incident_edge();
            const voronoi_edge<double>* start = edge;
            do {
                const voronoi_vertex<double>* v = edge->vertex0();
                if (v && (v->x() < xmin || v->x() > xmax || v->y() < ymin || v->y() > ymax)) {
                    is_external = true;
                    break;
                }
                edge = edge->next();
            } while (edge != start);
            
            if (!is_external) {
                double area = calculateVoronoiCellArea(*cell_it, vd);
                areas.push_back(area);
            }
        }
    }

    std::cout << "Areas: " << areas.size() << std::endl;
    std::cout << "Number of cells: " << vd.num_cells() << std::endl;

    // Calcula el área promedio
    double total_area = 0.0;
    for (double area : areas) {
        total_area += area;
    }
    double average_area = total_area / areas.size();

    // Establece el umbral de área a la mitad del área promedio
    double area_threshold = average_area * 0.8;

    double model_distance = sqrt(area_threshold/3.1416) * 2;

    //colorVoronoiCells(vd, area_threshold, xmin, xmax, ymin, ymax, model_distance, points);

    std::vector<double> normalized_areas;

    for (double area : areas) {
        double normalized_area = area / average_area;
        normalized_areas.push_back(normalized_area);
    }

    const double coefficient = std::pow(4, 4) / tgamma(4); // tgamma(n) es la función gamma en C++

    

    std::vector<double> dp_values;
    for (double normalized_area : normalized_areas) {
        double dp_value = dp(normalized_area, coefficient);
        dp_values.push_back(dp_value);
    }

    saveToCSV(normalized_areas, dp_values, "output.csv");


    std::vector<std::size_t> small_cells = findSmallCells(vd, area_threshold, xmin, xmax, ymin, ymax);
    visualizeVoronoiDiagram(vd, points);
    labelClusters(vd, points, model_distance);
    visualizeVoronoiDiagram(vd, points);
    
    connectNearbyCells(vd, points, model_distance, xmin,  xmax,  ymin, ymax);

    visualizeVoronoiDiagram(vd, points);
    auto clusters = collectClusters(vd, points);
    removeIsolatedSmallCells(vd, clusters, 8);
    
    
    visualizeVoronoiDiagram(vd, points);
    

    calculateClusterProperties(vd, points, "clusters.csv");

    return 0;
}





