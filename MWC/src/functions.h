#include <vector>
using namespace std;
// Function to search the minimum distance between a point that does not belong to the coalition and a array of distances and return the index of this
double minimum_distance(double distances[], bool congressman_set[], int n)
{
    double min = DBL_MAX, min_index;
    for (int v = 0; v < n; v++)
        if (congressman_set[v] == false && distances[v] < min)
            min = distances[v], min_index = v;
    return min_index;
}
// Function to evaluate the solutions and return the fitness value
double evaluate_solution(int *position, double **mat, int length)
{
    double sum = 0;
    for (size_t i = 0; i <= (length - 2); i++)
    {

        for (size_t j = i + 1; j <= (length - 1); j++)
        {
            sum = sum + mat[position[i]][position[j]];
        }
    }
    return sum;
}
// Function to create a initials solutions
void minimum_distance_edge(int *index_array, double *arrayDist, int n, int quorum)
{
    bool *congressman_set = new bool[n];
    for (size_t i = 0; i < n; i++)
    {
        congressman_set[i] = false;
    }
    for (size_t i = 0; i < quorum; i++)
    {
        int u = minimum_distance(arrayDist, congressman_set, n);
        congressman_set[u] = true;
        index_array[i] = u;
    }
}
// Function to calculate the distance between two points
double eucledian_distance(double x1, double y1, double x2, double y2)
{
    double calculation = pow(pow((x2 - x1), 2) + pow((y2 - y1), 2), 1 / (double)2);
    return calculation;
}
// Structure to store the information of the points
struct Point
{
    double x, y;
    int position;
    int index;
};
// Structure to store the distance respect to a centroid
struct Distance_vector
{
    double centroid_distance;
    int position;
};
// Structure to store the information of the points are possible to be selected to change
struct Possible_improvement
{
    double fitness;
    int index;
};
// Structure to store the information of the convex hull
struct Distance_hull
{
    double distance;
    int hull_index;
};
// Structure to store the initial solutions
struct Solutions
{
    double fitness;
    int *coalition_from_solution;
};
// Function to check the orientation
double orientation(Point p, Point q, Point r)
{
    double value = (q.y - p.y) * (r.x - q.x) -
                   (q.x - p.x) * (r.y - q.y);

    if (value == 0)
        return 0;               // collinear
    return (value > 0) ? 1 : 2; // clock or counterclock wise
}

// Prints convex hull of a set of n points.
// Source by GeeksforGeeks. (2022, August 5). Convex Hull | Set 1 (Jarvis’s Algorithm or Wrapping). Retrieved September 09, 2022, from https://www.geeksforgeeks.org/convex-hull-set-1-jarviss-algorithm-or-wrapping/
vector<Point> convexHull(Point points[], int n)
{
    // There must be at least 3 points
    // Initialize Result
    vector<Point> hull;
    if (n < 3)
        return hull;
    // Find the leftmost point
    int l = 0;
    for (int i = 1; i < n; i++)
        if (points[i].x < points[l].x)
            l = i;

    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    int p = l, q;
    do
    {
        // Add current point to result
        hull.push_back(points[p]);

        // Search for a point 'q' such that orientation(p, q,
        // x) is counterclockwise for all points 'x'. The idea
        // is to keep track of last visited most counterclock-
        // wise point in q. If any point 'i' is more counterclock-
        // wise than q, then update q.
        q = (p + 1) % n;
        for (int i = 0; i < n; i++)
        {
            // If i is more counterclockwise than current q, then
            // update q
            if (orientation(points[p], points[i], points[q]) == 2)
                q = i;
        }

        // Now q is the most counterclockwise with respect to p
        // Set p as q for next iteration, so that q is added to
        // result 'hull'
        p = q;

    } while (p != l); // While we don't come to first point
    // Print Result
    return hull;
}
// Boolean functions to check if a vector or array is sorted
bool vector_improvement_sort(Possible_improvement const &lvd, Possible_improvement const &rvd)
{
    return lvd.fitness < rvd.fitness;
}
bool vector_distance_hull_sort(Distance_hull const &lvd, Distance_hull const &rvd)
{
    return lvd.distance > rvd.distance;
}
bool vector_distance_sort(Distance_vector const &lvd, Distance_vector const &rvd)
{
    return lvd.centroid_distance < rvd.centroid_distance;
}
bool vector_initial_solutions_sort(Solutions const &lvd, Solutions const &rvd)
{
    return lvd.fitness < rvd.fitness;
}
bool array_sort(int const &lvd, int const &rvd)
{
    return lvd < rvd;
}
// Function to calculate the centroid of a set of points
void calculate_centroid(double *centroid, double **position_matrix, int *coalition, int quorum)
{
    centroid[0] = 0;
    centroid[1] = 0;
    for (size_t i = 0; i < quorum; i++)
    {
        centroid[0] = centroid[0] + position_matrix[coalition[i]][0];
        centroid[1] = centroid[1] + position_matrix[coalition[i]][1];
    }
    centroid[0] = centroid[0] / quorum;
    centroid[1] = centroid[1] / quorum;
}
// Funtion to calculate the distance between a point and a convex hull
vector<Distance_hull> distance_of_hull_to_point(vector<Point> hull, double *centroid)
{
    vector<Distance_hull> vector_distance_hull;
    for (size_t i = 0; i < hull.size(); i++)
    {
        vector_distance_hull.push_back(Distance_hull());
        vector_distance_hull[i].distance = eucledian_distance(hull[i].x, hull[i].y, centroid[0], centroid[1]);
        vector_distance_hull[i].hull_index = i;
    }
    sort(vector_distance_hull.begin(), vector_distance_hull.end(), &vector_distance_hull_sort);
    return vector_distance_hull;
}
// Function to calculate the distance
void distance_of_points_to_coalition(struct Distance_vector *distance_vector_minimum_winning_coalition, int *not_in_minimum_winning_coalition, int *coalition, double *centroid, double **position_matrix, int n, int quorum)
{
    for (int i = 0; i < (n - quorum); i++)
    {
        distance_vector_minimum_winning_coalition[i].position = not_in_minimum_winning_coalition[i];
        distance_vector_minimum_winning_coalition[i].centroid_distance = eucledian_distance(position_matrix[not_in_minimum_winning_coalition[i]][0], position_matrix[not_in_minimum_winning_coalition[i]][1], centroid[0], centroid[1]);
    }
    sort(distance_vector_minimum_winning_coalition, distance_vector_minimum_winning_coalition + (n - quorum), &vector_distance_sort);
}