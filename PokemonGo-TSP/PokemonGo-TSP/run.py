# coding:utf-8
from util import *
from args import *
from algorithm import *
from sympy.utilities.iterables import multiset_permutations

def Main(array, start, terminal):
    """
    main solving algorithm
    """
    n = array.size
    l = [start]
    for t in array:
        rand_theta = np.random.rand() * 2 * math.pi  # randomly initialization
        l.append(Point(Os[t].o.x + Os[t].r * math.cos(rand_theta),
                       Os[t].o.y + Os[t].r * math.sin(rand_theta)))
    l.append(terminal)
    rou = math.inf

    while abs(rou - Length(l)) > 10e-10:  # convergence
        i = 1
        while i < n + 1:
            # find the inverse map value, since the order is a map
            l[i] = GetoptPi(l[i - 1], l[i + 1], Os[array[i-1]])
            i = i + 1
        i = 2
        while i < n + 1:
            l[i] = GetoptPi(l[i - 1], l[i + 1], Os[array[i-1]])
            i = i + 2
        rou = Length(l)
    for x in range(len(l)):
        print(l[x].x,l[x].y)
    return rou, l


def main(scale, num):  # test with scale and num of circle
    # # Get argument
    args = parse_args()
    coordinates = scale * np.random.random(size=(num, 2))  # randomize a series of circle centers
    # Os initialization
    for x in range(coordinates.shape[0]):
        Os.append(Circle(Point(coordinates[x][0], coordinates[x][1]), radius))

    # Constant Definitions
    NUM_NEW_SOLUTION_METHODS = 3
    SWAP, REVERSE, TRANSPOSE = 0, 1, 2

    # Params Initial
    num_location = coordinates.shape[0]
    markov_step = args.markov_coefficient * num_location
    T_0, T, T_MIN = args.init_temperature, args.init_temperature, 1
    T_NUM_CYCLE = 1

    # States: New, Current and Best
    sol_new, sol_current, sol_best = (np.arange(num_location),) * 3
    cost_new, cost_current, cost_best = (float('inf'),) * 3

    # Record costs during the process
    costs = []

    # previous cost_best
    prev_cost_best = cost_best

    # counter for detecting how stable the cost_best currently is
    cost_best_counter = 0

    # Simulated Annealing
    while T > T_MIN and cost_best_counter < args.halt:
        for i in np.arange(markov_step):
            # Use three different methods to generate new solution
            # Swap, Reverse, and Transpose
            choice = np.random.randint(NUM_NEW_SOLUTION_METHODS)
            if choice == SWAP:
                sol_new = swap(sol_new)
            elif choice == REVERSE:
                sol_new = reverse(sol_new)
            elif choice == TRANSPOSE:
                sol_new = transpose(sol_new)
            else:
                print("ERROR: new solution method %d is not defined" % choice)

            # Get the total distance of new route and route itself
            cost_new, l_new = Main(sol_new, start, terminal)

            if accept(cost_new, cost_current, T):
                # Update sol_current
                sol_current = sol_new.copy()
                cost_current = cost_new

                if cost_new < cost_best:
                    l_best = l_new.copy()
                    sol_best = sol_new.copy()
                    cost_best = cost_new
            else:
                sol_new = sol_current.copy()

        # Lower the temperature
        alpha = 1 + math.log(1 + T_NUM_CYCLE)
        T = T_0 / alpha
        costs.append(cost_best)

        # Increment T_NUM_CYCLE
        T_NUM_CYCLE += 1

        # Detect stability of cost_best
        if isclose(cost_best, prev_cost_best, abs_tol=1e-12):
            cost_best_counter += 1
        else:
            # Not stable yet, reset
            cost_best_counter = 0

        # Update prev_cost_best
        prev_cost_best = cost_best

        # Monitor the temperature & cost
        print("Temperature:", "%f°C" % T,
              " Distance:", "%f" % cost_best,
              " Optimization Threshold:", "%d" % cost_best_counter)
    # Show final cost & route
    print("Final Distance:", costs[-1])
    print("Best Route", sol_best)

    # Plot cost function and TSP-route
    # #TODO test it
    # sol_best=np.array([2,3,1,0])
    # _,l_best=Main(sol_best,start,terminal)
    plot(sol_best, l_best, coordinates)



def plot(path, lk, points, costs=None):
    '''
    draw the plots
    '''

    # Change figure size
    plt.figure(figsize=(15, 6))

    '''
    Plot Cost Function
    '''
    if costs:
        plt.subplot(121)
        curve, = plt.plot(np.array(costs), label='Distance(m)')
        plt.ylabel("Distance")
        plt.xlabel("Iteration")
        plt.grid(True)
        plt.legend()
        cost = str("%.2f" % round(costs[-1], 2))
        plt.title("Final Distance: " + cost)

    '''
    Plot TSP Route
    '''
    figure, axes = plt.subplots()
    points = points.tolist()
    x = []
    y = []
    for _k in lk:
        x.append(_k.x)
        y.append(_k.y)
    for i in Os:
        Drawing_uncolored_circle = plt.Circle((i.o.x, i.o.y),
                                              i.r,
                                              fill=False)

        axes.set_aspect(1)
        axes.add_artist(Drawing_uncolored_circle)
    for j in range(len(x)):
        print(x[j],y[j])
    plt.plot(x, y, 'c-', label='Route')
    plt.plot(x, y, 'bo', label='Location')

    # Avoid scientific notation
    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    # Set axis too slightly larger than the set of x and y
    plt.xlim(min(x) * 0.99999, max(x) * 2.00001)
    plt.ylim(min(y) * 0.99999, max(y) * 2.00001)
    plt.xlabel("axis0")
    plt.ylabel("axis1")
    plt.title("TSP Route Visualization")
    plt.grid(True)
    plt.show()


"""global variable"""
radius = 200
Os = []
start = Point(0, 0)
terminal = Point(0, 0)
"""
for test(II)
"""
# Os = [Circle(Point(8,8),radius),Circle(Point(4,8),radius),Circle(Point(12,8),radius),
#       Circle(Point(6,8+2*math.sqrt(3)),radius),Circle(Point(6,8-2*math.sqrt(3)),radius),
#       Circle(Point(10,8+2*math.sqrt(3)),radius),Circle(Point(10,8-2*math.sqrt(3)),radius)]


if __name__ == "__main__":
    # init parameters
    scale = 1000
    num = 12
    main(scale, num)
