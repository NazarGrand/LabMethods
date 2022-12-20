using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp1
{
    class Ortogonalization
    {
        private double[,] A = new double[3, 3] { { 2.16, -3.18, 1.26 }, { -3.18, 0.63, -2.73 }, {1.26, -2.73, 3.15 } };
        private double[] b = new double[3] { 1.83, 0.54, 1.72 };

        private double[] x1 = new double[3] { 1, 0, 0 }, x2 = new double[3], x3 = new double[3];
        private double a12, a13, a23;

        private double[] e1 = new double[3] { 1, 0, 0 }, e2 = new double[3] { 0, 1, 0 }, e3 = new double[3] { 0, 0, 1 };
        private double C1, C2, C3;

        public double[] X = new double[3];
        public double[] q = new double[3];

        public double DeterminantMatrix()
        {
            return A[0, 0] * A[1, 1] * A[2, 2] + A[0, 1] * A[1, 2] * A[2, 0] + A[1, 0] * A[2, 1] * A[0, 2] -
                   A[0, 2] * A[1, 1] * A[2, 0] - A[1, 0] * A[0, 1] * A[2, 2] - A[2, 1] * A[1, 2] * A[0, 0];
        }

        public double[] MatrixMultVector(double[] vec)
        {
            double[] res = new double[3];

            for (int i = 0; i < A.GetLength(0); i++)
                for (int j = 0; j < A.GetLength(1); j++)
                    res[i] += A[i, j] * vec[j];

            return res;
        }

        public double ScalarDobutoc(double[] vec1, double[] vec2)
        {
            return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
        }

        public double[] MultVectorNumber(double[] vec, double num)
        {
            return new double[3] { vec[0] * num, vec[1] * num, vec[2] * num };
        }

        public double[] AdditionVectors(double[] vec1, double[] vec2)
        {
            return new double[3] { vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2] };
        }

        public void FindSolution()
        {
            if (DeterminantMatrix() != 0)
            {
                double[] vec1 = MatrixMultVector(x1);
                a12 = -ScalarDobutoc(vec1, e2) / ScalarDobutoc(vec1, x1);
                x2 = AdditionVectors(e2, MultVectorNumber(x1, a12));

                a13 = -ScalarDobutoc(vec1, e3) / ScalarDobutoc(vec1, x1);

                double[] vec2 = MatrixMultVector(x2);
                a23 = -ScalarDobutoc(vec2, e3) / ScalarDobutoc(vec2, x2);
                x3 = AdditionVectors(AdditionVectors(e3, MultVectorNumber(x1, a13)), MultVectorNumber(x2, a23));

                C1 = ScalarDobutoc(b, x1) / ScalarDobutoc(MatrixMultVector(x1), x1);
                C2 = ScalarDobutoc(b, x2) / ScalarDobutoc(MatrixMultVector(x2), x2);
                C3 = ScalarDobutoc(b, x3) / ScalarDobutoc(MatrixMultVector(x3), x3);

                X = AdditionVectors(AdditionVectors(MultVectorNumber(x1, C1), MultVectorNumber(x2, C2)), MultVectorNumber(x3, C3));
                Inconsistencies();
            }
            else
                Console.WriteLine("Determinant = 0");
        }

        public void Inconsistencies()
        {
            q[0] = X[0] * A[0, 0] + X[1] * A[0, 1] + X[2] * A[0, 2] - b[0];
            q[1] = X[0] * A[1, 0] + X[1] * A[1, 1] + X[2] * A[1, 2] - b[1];
            q[2] = X[0] * A[2, 0] + X[1] * A[2, 1] + X[2] * A[2, 2] - b[2];
        }
    }

    class EidelMethod
    {
        private double[,] A = new double[3, 3] { { 2.16, -3.18, 1.26 }, { -3.18, 0.63, -2.73 }, { 1.26, -2.73, 3.15 } };
        private double[] b = new double[3] { 1.83, 0.54, 1.72 };
        private double[,] newA = new double[3, 3];
        private double[] newB = new double[3];

        public double[] x = new double[3];
        public int numberIteration = 0;
        public double[] q = new double[3];

        public void MatrixTransformation()
        {
           for(int i = 0; i < A.GetLength(1); i++)
            {
                newA[0, i] = 2 * A[1, i] + A[2, i]; //2*II + III
                newA[1, i] = A[0, i] + A[1, i]; //I + II
                newA[2, i] = A[2, i] - A[0, i]; //III - I
            }

            newB[0] = 2 * b[1] + b[2];
            newB[1] = b[0] + b[1];
            newB[2] = b[2] - b[0];
        }

        public void SolutionMethod()
        {
            MatrixTransformation();

            x[0] = newB[0]; x[1] = newB[1]; x[2] = newB[2];

            double[] tempX = new double[3];

            do
            {
                tempX[0] = x[0];
                tempX[1] = x[1];
                tempX[2] = x[2];
                x[0] = (newB[0] - newA[0, 1] * x[1] - newA[0, 2] * x[2]) / newA[0, 0];
                x[1] = (newB[1] - newA[1, 0] * x[0] - newA[1, 2] * x[2]) / newA[1, 1];
                x[2] = (newB[2] - newA[2, 0] * x[0] - newA[2, 1] * x[1]) / newA[2, 2];
                numberIteration++;
            }
            while ((Math.Abs(x[0] - tempX[0]) > 0.001) && (Math.Abs(x[1] - tempX[1]) > 0.001) &&
                   (Math.Abs(x[2] - tempX[2]) > 0.001));

            Inconsistencies();
        }

        public void Inconsistencies()
        {
            q[0] = x[0] * newA[0, 0] + x[1] * newA[0, 1] + x[2] * newA[0, 2] - newB[0];
            q[1] = x[0] * newA[1, 0] + x[1] * newA[1, 1] + x[2] * newA[1, 2] - newB[1];
            q[2] = x[0] * newA[2, 0] + x[1] * newA[2, 1] + x[2] * newA[2, 2] - newB[2];
        }
    }

    internal class Program
    {
        static void Main(string[] args)
        {
            Ortogonalization ortogonalization = new Ortogonalization();
            ortogonalization.FindSolution();
            Console.WriteLine("Orthogonalization method:");
            Console.WriteLine($"x1 = {ortogonalization.X[0]}");
            Console.WriteLine($"x2 = {ortogonalization.X[1]}");
            Console.WriteLine($"x3 = {ortogonalization.X[2]}");
            Console.WriteLine("Inconsistencies:");
            Console.WriteLine($"q1 = {ortogonalization.q[0]}");
            Console.WriteLine($"q2 = {ortogonalization.q[1]}");
            Console.WriteLine($"q3 = {ortogonalization.q[2]}");

            EidelMethod eidelMethod = new EidelMethod();
            eidelMethod.SolutionMethod();
            Console.WriteLine("\n\nEidel method:");
            Console.WriteLine($"x1 = {eidelMethod.x[0]}");
            Console.WriteLine($"x2 = {eidelMethod.x[1]}");
            Console.WriteLine($"x3 = {eidelMethod.x[2]}");
            
            Console.WriteLine("Inconsistencies:");
            Console.WriteLine($"q1 = {eidelMethod.q[0]}");
            Console.WriteLine($"q2 = {eidelMethod.q[1]}");
            Console.WriteLine($"q3 = {eidelMethod.q[2]}");
            Console.WriteLine($"Number of iteration = {eidelMethod.numberIteration}");

            Console.ReadKey();
        }
    }
}
