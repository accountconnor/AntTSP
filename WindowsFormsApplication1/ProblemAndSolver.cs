using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;


namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            string[] results = new string[3];

            /////////MY CODE BELOW/////////////////////////////////////////////////////////////////////////////

            Stopwatch timer = new Stopwatch();
            int count = 0;
            int pruneCount = 0;
            int bssfUpdateCount = 0;
            int greatestTotalStates = 0;
            timer.Start();

            List<List<double>> graph = generateGraph();
            Tuple<double,List<int>> bestPath = findUpperBound(graph);
            List<State> stateList = new List<State>();
            stateList.Add(new State(graph));
            State minLowerBoundState = new State();
            List<State> newStates = new List<State>();
            int minLowerBoundIndex = 0;

            
            /* This next section of code runs once for every State that's generated but not pruned.
                This totals to
            */
            while (stateList.Count > 0 && timer.Elapsed.TotalMilliseconds < time_limit)
            {
                if(stateList.Count > greatestTotalStates)
                {
                    greatestTotalStates = stateList.Count;
                }
                minLowerBoundIndex = findMinLowerBoundIndex(stateList);
                minLowerBoundState = stateList[minLowerBoundIndex];
                newStates = expandState(minLowerBoundState);
                count += newStates.Count;
                stateList.Remove(minLowerBoundState);
                for (int i = 0; i < newStates.Count; i++)
                {
                    switch (newStates[i].finished(bestPath.Item1))
                    {
                        case 1:
                            double candidate = newStates[i].getCurrentCost() + newStates[i].getGraph()[newStates[i].getCurrentNodeIndex()][0];
                            if (candidate < bestPath.Item1) {
                                bssfUpdateCount++;
                                bestPath = Tuple.Create<double, List<int>>(candidate, newStates[i].getRoute());
                            }
                            break;
                        case 0:
                            stateList.Add(newStates[i]);
                            break;
                        case 2:
                            pruneCount++;
                            break;
                    }
                }
            }

            ArrayList foundBSSF = new ArrayList();
            for(int i = 0; i < bestPath.Item2.Count; i++){
                foundBSSF.Add(Cities[bestPath.Item2[i]]);
            }
            bssf = new TSPSolution(foundBSSF);

            timer.Stop();

            results[COST] = bestPath.Item1.ToString();
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            Console.WriteLine(greatestTotalStates);
            Console.WriteLine(bssfUpdateCount);
            Console.WriteLine(pruneCount);

            return results;
        }

        public class State
        {
            private double currentCost = 0;
            private double lowerBound = 0;
            private int currentNodeIndex = 0;
            private List<int> route = new List<int>();
            private List<List<double>> graph = new List<List<double>>();
            
            public State() {
                route.Add(0);
            }

            public State(bool addZero)
            {
                if (addZero)
                    route.Add(0);
            }

            public State(List<List<double>> nodeGraph)
            {
                graph = nodeGraph;
                route.Add(0);
            }

            public double getLowerBound()
            {
                return lowerBound;
            }
      
            // This is so I can make a deep copy, rather than a shallow one.
            public State clone()
            {
                State result = new State(false);

                result.currentCost = currentCost;
                result.lowerBound = lowerBound;

                for(int i = 0; i < route.Count; i++)
                {
                    result.addToRoute(route[i]);
                }
                result.currentNodeIndex = currentNodeIndex;
                for (int i = 0; i < graph.Count; i++)
                {
                    result.graph.Add(new List<double>());
                    for (int j = 0; j < graph.Count; j++)
                    {
                        result.graph[i].Add(graph[i][j]);
                    }
                }

                return result;
            }

            public double getCurrentCost()
            {
                return currentCost;
            }

            // O(n^2) where n is the length of the graph.
            public void updateLowerBound()
            {
                double result = currentCost;
                int index = 0;
                for(int i = 0; i < graph.Count; i++)
                {
                    if (!route.Contains(i))
                    {
                        index = 0;
                        for (int j = 1; j < graph[i].Count && i != j; j++)
                        {
                            if (graph[i][j] < graph[i][index])
                            {
                                index = j;
                            }
                        }
                        if (index != 0)
                        {
                            result += graph[i][index] - 1;
                        }
                        else
                        {
                            result += (int)(graph[i][index] / 2);
                        }
                    }
                }
                lowerBound = result;
            }

            public void addToCurrentCost(double toAdd)
            {
                currentCost += toAdd;
            }

            public List<int> getRoute()
            {
                return route;
            }

            public List<List<double>> getGraph()
            {
                return graph;
            }

            public void addToRoute(int input)
            {
                route.Add(input);
            }
            
            /*
               0) Not finished
               There are two ways to be finished:
               1) You have visited all the nodes, which means the route is as big as there are nodes.
               2) You reach a node with no viable outnodes or the lowerbound is greater than the upper bound.
               If it's for the first reason, then you need to compare the route and cost with the upperBound
               and see if it needs to be replaced. If it's for the second reason, then the state must be pruned.
            */
            public int finished(double upperBound)
            {
                if(lowerBound > upperBound)
                {
                    return 2;
                }
                if(route.Count >= graph.Count)
                {
                    return 1;
                }
                int result = 2;
                for(int i = 1; i < graph.Count; i++)
                {
                    if (!double.IsInfinity(graph[currentNodeIndex][i]))
                    {
                        result = 0;
                    }
                }
                return result;
            }

            public void clearRowAndColumn(int row, int col)
            {
                for(int i = 0; i < graph.Count; i++)
                {
                    graph[row][i] = double.PositiveInfinity;
                    graph[i][col] = double.PositiveInfinity;
                }

            }

            public int getCurrentNodeIndex()
            {
                return currentNodeIndex;
            }

            public void setCurrentNodeIndex(int newIndex)
            {
                currentNodeIndex = newIndex;
            }
        }

        // O(n^3) -> where n is the number of unvisited nodes.
        // The first n comes because the number of expansions is
        // equal to the number of unvisited nodes that you can reach
        // from our input state.
        // The second two n's come together because of the
        // updateLowerBound() method.
        public List<State> expandState(State input)
        {
            List<State> result = new List<State>();
            int current = input.getCurrentNodeIndex();
            List<double> node = input.getGraph()[current];
            State addableState = new State();
            for (int i = 1; i < node.Count; i++)
            {
                if (!double.IsPositiveInfinity(node[i])) {
                    addableState = input.clone();
                    addableState.addToCurrentCost(node[i]);
                    addableState.setCurrentNodeIndex(i);
                    addableState.clearRowAndColumn(current, i);
                    addableState.updateLowerBound();
                    addableState.addToRoute(i);
                    result.Add(addableState);
                }
            }
            return result;
        }

        // O(n) where n is the number of states.
        public int findMinLowerBoundIndex(List<State> input)
        {
            int result = 0;
            for(int i = 0; i < input.Count; i++)
            {
                if(input[i].getLowerBound() < input[result].getLowerBound())
                {
                    result = i;
                }
            }
            return result;
        }

        public List<List<double>> generateGraph()
        {
            List<List<double>> result = new List<List<double>>();
            for (int i = 0; i < Cities.Length; i++)
            {
                result.Add(new List<double>());
                for(int j = 0; j < Cities.Length; j++)
                {
                    result[i].Add(i!=j?Cities[i].costToGetTo(Cities[j]):double.PositiveInfinity);
                }
            }
            return result;
        }

        public Tuple<double, List<int>> findUpperBound(List<List<double>> graph)
        {
            Tuple<double, List<int>> resultPath = new Tuple<double, List<int>>(double.PositiveInfinity,new List<int>());
            Tuple<double, List<int>> tempPath = new Tuple<double, List<int>>(double.PositiveInfinity, new List<int>());
            for (int i = 0; i < 10; i++)
            {
                tempPath = findRandomPath(graph);
                if(tempPath != null)
                {
                    if(tempPath.Item1 < resultPath.Item1)
                    {
                        resultPath = tempPath;
                    }
                }
            }
            return resultPath;
        }

        // Generate random tours through the graph to get an idea as to what
        // the upper bound might be. It's O(n) where n is the size of the graph,
        // but there's a big constant of 100000 attached to that.
        public Tuple<double, List<int>> findRandomPath(List<List<double>> graph)
        {
            List<List<Boolean>> paths = new List<List<Boolean>>();
            for (int i = 0; i < graph.Count; i++)
            {
                paths.Add(new List<Boolean>());
                for (int j = 0; j < graph[i].Count; j++)
                {
                    paths[i].Add(graph[i][j] != double.PositiveInfinity);
                }
            }

            List<int> result = new List<int>();
            for (int i = 0; i < graph.Count; i++)
            {
                result.Add(i);
            }
            bool done = false;
            int count = 100000;
            while (!done && count > 0)
            {
                done = true;
                count--;
                result = shuffle(result);
                for (int i = 0; i < paths.Count; i++)
                {
                    if (!paths[i][result[i]])
                    {
                        done = false;
                        i = paths.Count;
                    }
                }
            }
            if (count == 0)
            {
                return null;
            }
            double cost = 0;
            for (int i = 0; i < graph.Count; i++)
            {
                cost += graph[i][result[i]];
            }

            return Tuple.Create<double, List<int>>(cost, result);
        }

        public List<int> shuffle(List<int> input)
        {
            List<int> result = new List<int>();
            Random rnd = new Random();
            while (input.Count > 0)
            {
                int index = rnd.Next(0, input.Count);
                result.Add(input[index]);
                input.RemoveAt(index);
            }
            return result;
        }
        /*
        public Tuple<double,List<List<double>>> residualize(List<List<double>> input)
        {
            double resultTotal = 0;
            List<double> toDo = new List<double>();

            for(int i = 0; i < input.Count; i++)
            {
                toDo.Add(i);
            }

            Tuple<int, double> listBest = new Tuple<int, double>(0,0);
            for(int i = 0; i < input.Count; i++)
            {
                listBest = listMin(input[i]);
                for(int j = 0; j < input[i].Count; j++)
                {
                    input[i][j] -= listBest.Item2;
                }
                toDo.Remove(listBest.Item1);
                resultTotal += listBest.Item2;
            }
            for (int i = 0; i < toDo.Count; i++)
            {

                listBest = listMin(getColumn(input,i));
                for (int j = 0; j < input[i].Count; j++)
                {
                    input[j][i] -= listBest.Item2;
                }
                resultTotal += listBest.Item2;

            }
            return Tuple.Create(resultTotal, input);
        }

        List<double> getColumn(List<List<double>> matrix, int index)
        {
            List<double> result = new List<double>();
            for(int i = 0; i < matrix.Count; i++)
            {
                result.Add(matrix[i][index]);
            }
            return result;
        }

        // returns the index and value of the min.
        Tuple<int, double> listMin(List<double> input)
        {
            double resultValue = 0;
            int resultIndex = 0;
            for(int i = 0; i < input.Count; i++)
            {
                if (input[i] < resultValue)
                {
                    resultValue = input[i];
                    resultIndex = i;
                }
            }
            return Tuple.Create<int, double>(resultIndex, resultValue);
        }
        */

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];
            
            // TODO: Add your implementation for a greedy solver here.

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }

        public string[] fancySolveProblem()
        {
            string[] results = new string[3];


            Stopwatch timer = new Stopwatch();
            int count = 0;
            int bssfUpdateCount = 0;
            int greatestTotalStates = 0;
            double pheremoneRate = 1.0;
            int iterations = 10000;
            SendAntResult sendAntResult;
            timer.Start();
            
            List<List<double>> graph = generateGraph();

            List<List<double>> pheremones = new List<List<double>>();
            for (int i = 0; i < graph.Count; i++)
            {
                pheremones.Add(new List<double>());
                for(int j = 0; j < graph.Count; j++)
                {
                    pheremones[i].Add(0.0);
                }
            }

            for (int i = 0; i < iterations; i++)
            {
                sendAntResult = sendAnt(graph, pheremones, pheremoneRate);
                pheremones = sendAntResult.NewPheremones;
                //TODO: update bssf if new cost is less than old cost.
            }
            
            //TODO: find result route from pheremone trail. Or, we could just send the bssf.


                ///////////////////////////////////////////////////////////////

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }

        // Generate random tours through the graph to get an idea as to what
        // the upper bound might be. It's O(n) where n is the size of the graph,
        // but there's a big constant of 100000 attached to that.
        public SendAntResult sendAnt(List<List<double>> graph, List<List<double>> pheremones, double pheremoneRate)
        {
            List<int> resultRoute = new List<int>();
            double cost = 0.0;

            int currentNode = 0;
            int nextEdge = 0;
            while(nextEdge == 0){
                nextEdge = selectEdge(pheremones[currentNode], graph[currentNode]);
                if (nextEdge == -1)
                {
                    return new SendAntResult(double.PositiveInfinity, null, pheremones);
                }

                cost += graph[currentNode][nextEdge];
                graph = clearRowAndColumn(graph, currentNode, nextEdge);
                currentNode = nextEdge;
            }

            
            //TODO: modify pheremones
           

            SendAntResult result = new SendAntResult(cost, resultRoute, pheremones);

            return result;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pheremones"></param>
        /// <param name="edges"></param>
        /// <returns></returns>
        int selectEdge(List<double> pheremones, List<double> edges)
        {
            int result = 0;

            double total = 0;
            double newTotal = 0;
            double probabilityTotal = 0;
            Dictionary<double, double> probability = new Dictionary<double, double>();
            Random rnd = new Random();

            foreach(double edge in edges)
            {
                if (edge != Double.PositiveInfinity)
                {
                    total += edge;
                }
            }

            foreach (double edge in edges)
            {
                if (edge != Double.PositiveInfinity)
                {
                    double num = total - edge;
                    newTotal += num;
                }
            }

            foreach(double edge in edges)
            {
                if(edge != Double.PositiveInfinity)
                {
                    double num = total - edge;
                    double prob = (num / newTotal) * 100;
                    probabilityTotal += prob;
                    probability.Add(edge, probabilityTotal);
                }
            }
            double random = rnd.NextDouble() * (100);

            double returnEdge = 0;
            foreach(var prob in probability)
            {
                if(random < prob.Value)
                {
                    returnEdge = prob.Key;
                    break;
                }
            }

            result = edges.IndexOf(returnEdge);
            return result;
        }

        public class SendAntResult
        {
            public double Cost { get; set; }
            public List<int> Route { get; set; }
            public List<List<double>> NewPheremones { get; set; }
            public SendAntResult(double cost, List<int> route, List<List<double>> newPheremones) {
                Cost = cost;
                Route = route;
                NewPheremones = newPheremones;
            }
        }

        public List<List<double>> clearRowAndColumn(List<List<double>> input, int row, int col)
        {
            for (int i = 0; i < input.Count; i++)
            {
                input[row][i] = double.PositiveInfinity;
                input[i][col] = double.PositiveInfinity;
            }
            return input;
        }

        #endregion
    }

}
