#include<iostream>
#include<fstream>
#include<unordered_set>
#include<utility>
#include<vector>
#include<random>
#include<cassert>
#include<algorithm>
#include<numeric>

using namespace std;

double** create_matrix(size_t a, size_t b)
{
	double ** m = new double *[a];
	m[0] = new double[a * b];
	for (size_t i = 1; i != a; ++i)
		m[i] = m[i - 1] + b;
	return m;
}


void free_matrix(double** m, size_t a, size_t b)
{
	delete[] m[0];
	delete[] m;
}

void read_matrix(std::istream& input, double ** m, size_t a, size_t b)
{
	for (size_t i = 0; i < a; ++i)
		for (size_t j = 0; j < b; ++j)
			input >> m[i][j];
}

void print_matrix(std::ostream& output, double ** m, size_t a, size_t b)
{
	for (size_t i = 0; i < a; ++i)
	{
		for (size_t j = 0; j < b; ++j)
			output << m[i][j] << ' ';
		output << endl;
	}
}

void clear_matrix(double ** m, size_t a, size_t b)
{
	memset(m[0], 0, a * b * sizeof(double));
}

vector<vector<size_t>> find_cycles(double ** m, size_t a)
{
	return{ { 0, 1, 2, 3 }, { 0, 1, 2, 3, 4} };
}

double d_max(double ** m, size_t a, size_t b)
{
	double dmax = 0;
	for (size_t i = 0; i < a; ++i)
	{
		double deg = 0;
		for (size_t j = 0; j < a; ++j)
		{
			deg += m[i][j];
		}
		dmax = max(dmax, deg);
	}
	return dmax;
}

void multiply_scalar(double ** m, size_t a, size_t b, double scalar)
{
	for (size_t i = 0; i < a; ++i)
		for (size_t j = 0; j < b; ++j)
			m[i][j] *= scalar;
}

void apply_laplacian(double ** m, size_t a)
{
	double** d = create_matrix(a,a);
	for (size_t i = 0; i < a; ++i)
	{
		double deg = 0;
		for (size_t j = 0; j < a; ++j)
		{
			deg += m[i][j];
			if (i != j)
				d[i][j] = 0;
		}
		d[i][i] = deg;
	}
	multiply_scalar(m, a, a, -1);
	for (size_t i = 0; i < a; ++i)
		for (size_t j = 0; j < a; ++j)
			m[i][j] += d[i][j];
}

void add_scalar(double ** m, size_t a, size_t b, double scalar)
{
	for (size_t i = 0; i < a; ++i)
		for (size_t j = 0; j < b; ++j)
			m[i][j] += scalar;
}


void vector_multiply(double ** B, size_t a, size_t b, vector<double>& x)
{
	vector<double> new_x(x.size());

	for (size_t i = 0; i < new_x.size(); ++i)
		for (size_t j = 0; j < a; ++j)
			new_x[i] += B[i][j] * x[j];

	x.assign(new_x.cbegin(),new_x.cend());
}

double calc_amplitude(const vector<double>& x1)
{
	double res = 0;
	for (size_t i = 0; i < x1.size(); ++i)
	{
		double curr = fabs(x1[i]);
		res = max(res, curr);
	}
	return res;
}


double calc_deviations(double ** B, size_t a, const vector<double>& x0)
{
	//double avg = accumulate(x0.cbegin(), x0.cend(), 0) / x0.size();

	vector<double> x1(x0);

	double max_amp = 0;

	for (size_t i = 0; i < 100; ++i)
	{
		vector_multiply(B, a, a, x1);

		double amp = calc_amplitude(x1);
		//cout << amp << endl;
		max_amp = max(max_amp, amp);
	}

	return max_amp;
}




bool matrix_multiply(double**& m_first, size_t& rows_first, size_t& cols_first, double const * const * m_second, size_t rows_second, size_t cols_second)
{
	if (!(cols_first == rows_second))
		return false;

	double** res = create_matrix(rows_first, cols_second);
	memset(res[0], 0, rows_first * cols_second * sizeof(double));

	for (size_t i = 0; i < rows_first; ++i)
		for (size_t j = 0; j < cols_second; ++j)
			for (size_t k = 0; k < cols_first; ++k)
				res[i][j] += m_first[i][k] * m_second[k][j];


	free_matrix(m_first, rows_first, cols_first);
	m_first = res;
	cols_first = cols_second;

	return true;
}

void gen_control_matrix(double** B, size_t a, size_t b)
{
	double at = 1.0 / (d_max(B,a,a) + 1) / 10;
	//double at = 1.0 / 3;

	multiply_scalar(B, a, b, at);
	apply_laplacian(B, a);
	multiply_scalar(B, a, b, -1);
}

int main(int argc, char *argv[])
{
	const size_t ITERATIONS = 10;
	const size_t CYCLES = 10;

	if (argc != 2)
	{
		cout << "Error! Wrong number of arguments." << endl;
		return 1;
	}

	minstd_rand rng;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	ifstream input(argv[1]);
	size_t rows_current = 0;
	input >> rows_current;

	assert(rows_current > 0);

	vector<double> x0(rows_current);
	double max_x0 = dis(gen) * 5000;
	for (size_t i = 0; i < x0.size(); ++i)
		x0[i] = -max_x0 + 2 * (max_x0 / x0.size() * (i + 1));

	double** A = create_matrix(rows_current, rows_current);
	read_matrix(input, A, rows_current, rows_current);
	print_matrix(cout, A, rows_current, rows_current);
	vector<vector<size_t>> cycles = find_cycles(A, rows_current);

	double** B = create_matrix(rows_current, rows_current);
	double** B_best = create_matrix(rows_current, rows_current);
	double dev_best = 1.79769e+308;

	for (size_t iter = 0; iter < ITERATIONS; ++iter)
	{
		clear_matrix(B, rows_current, rows_current);

		for (size_t c = 0; c < CYCLES; ++c)
		{
			size_t curr_cycle = rng() % cycles.size();

			double delta = dis(gen) * 10;
			for (size_t i = 0; i < cycles[curr_cycle].size() - 1; ++i)
			{
				B[cycles[curr_cycle][i]][cycles[curr_cycle][i + 1]] += delta;
			}
			B[cycles[curr_cycle][cycles[curr_cycle].size() - 1]][cycles[curr_cycle][0]] += delta;

		}

		gen_control_matrix(B, rows_current, rows_current);

		double deviations = calc_deviations(B, rows_current, x0);
		cout << deviations << endl;

		if (deviations < dev_best)
		{
			for (size_t i = 0; i < rows_current; ++i)
				for (size_t j = 0; j < rows_current; ++j)
					B_best[i][j] = B[i][j];
			dev_best = deviations;
		}

	}

	ofstream output("res.txt");
	print_matrix(output, B_best, rows_current, rows_current);

	system("PAUSE");

	return 0;
}