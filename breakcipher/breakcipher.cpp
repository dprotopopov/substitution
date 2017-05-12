// Атака на шифр простой замены.
// Дмитрий Протопопов (protopopov@narod.ru)

#include "stdafx.h"
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace std;

namespace substitution
{
	// Класс частотной статистики триграмм
	class trigram
	{
	public:
		string alphabet; // Символы алфавита
		vector<long> data; // Частотная статистика триграмм
		int count; // Количество символов в алфавите
		long length; // Размер репрезентативной выборки

		explicit trigram(const string &fileName)
		{
			ifstream file(fileName);
			for (auto ch = file.get(); ch != EOF; ch = file.get())
				if (alphabet.find(ch) == string::npos)
					alphabet.append(1, ch);
			count = alphabet.length();
			data.resize(count*count*count);
			file.clear();
			file.seekg(0);
			auto index0 = 0;
			auto index1 = 0;
			for (auto ch = file.get(); ch != EOF; ch = file.get())
			{
				auto index2 = alphabet.find(ch);
				if (++length >= 3)
				{
					data[index0 + index1 * count + index2 * count * count]++;
					index0 = index1;
					index1 = index2;
				}
			}
			file.close();
		}
	};

	// Класс частотной статистики биграмм
	class bigram
	{
	public:
		string alphabet; // Символы алфавита
		vector<long> data; // Частотная статистика биграмм
		int count; // Количество символов в алфавите
		long length; // Размер репрезентативной выборки

		explicit bigram(const string &fileName)
		{
			ifstream file(fileName);
			for (auto ch = file.get(); ch != EOF; ch = file.get())
				if (alphabet.find(ch) == string::npos)
					alphabet.append(1, ch);
			count = alphabet.length();
			data.resize(count*count);
			file.clear();
			file.seekg(0);
			auto index0 = 0;
			for (auto ch = file.get(); ch != EOF; ch = file.get())
			{
				auto index1 = alphabet.find(ch);
				if (++length >= 2)
				{
					data[index0 + index1 * count]++;
					index0 = index1;
				}
			}
			file.close();
		}
		explicit bigram(const trigram &t)
		{
			alphabet = t.alphabet;
			length = t.length;
			count = t.count;
			data.resize(count*count);
			for (auto i = 0; i < count*count; i++)
			{
				data[i] = t.data[i];
				for (auto j = 1; j < count; j++)
				{
					data[i] += t.data[i + j*count*count];
				}
			}
		}
	};

	// Класс частотной статистики символов
	class unogram
	{
	public:
		string alphabet; // Символы алфавита
		vector<long> data; // Частотная статистика символов
		int count; // Количество символов в алфавите
		long length; // Размер репрезентативной выборки

		explicit unogram(const string &fileName)
		{
			ifstream file(fileName);
			for (auto ch = file.get(); ch != EOF; ch = file.get())
				if (alphabet.find(ch) == string::npos)
					alphabet.append(1, ch);
			count = alphabet.length();
			data.resize(count);
			file.clear();
			file.seekg(0);
			for (auto ch = file.get(); ch != EOF; ch = file.get())
			{
				length++;
				data[ch]++;
			}
			file.close();
		}
		explicit unogram(const bigram &b)
		{
			alphabet = b.alphabet;
			length = b.length;
			count = b.count;
			data.resize(count);
			for (auto i = 0; i < count; i++)
			{
				data[i] = b.data[i];
				for (auto j = 1; j < count; j++)
				{
					data[i] += b.data[i + j*count];
				}
			}
		}
	};

	// Расчёт корреляции по частотной статистике символов
	double fitness(const vector<int> &p1, const vector<int> &p2, const unogram &u1, const unogram &u2)
	{
		double s = 0;
		auto count = min(p1.size(), p2.size());
		for (auto i = 0; i < count; i++)
		{
			double x = (p1[i] < u1.count) ? u1.data[p1[i]] : 0;
			double y = (p2[i] < u2.count) ? u2.data[p2[i]] : 0;
			s += x*y;
		}
		return s;
	}

	// Расчёт корреляции по частотной статистике биграмм
	double fitness(const vector<int> &p1, const vector<int> &p2, const bigram &b1, const bigram &b2)
	{
		double s = 0;
		auto count = min(p1.size(), p2.size());
		auto count1 = b1.count;
		auto count2 = b2.count;
		for (auto i = 0; i < count; i++)
		{
			for (auto j = 0; j < count; j++)
			{
				double x = (p1[i] < count1&&p1[j] < count1) ? b1.data[p1[i] + p1[j] * count1] : 0;
				double y = (p2[i] < count2&&p2[j] < count2) ? b2.data[p2[i] + p2[j] * count2] : 0;
				s += x*y;
			}
		}
		return s;
	}

	// Расчёт корреляции по частотной статистике триграмм
	double fitness(const vector<int> &p1, const vector<int> &p2, const trigram &t1, const trigram &t2)
	{
		double s = 0;
		auto count = min(p1.size(), p2.size());
		auto count1 = t1.count;
		auto count2 = t2.count;
		auto count12 = count1*count1;
		auto count22 = count2*count2;

		for (auto i = 0; i < count; i++)
		{
			for (auto j = 0; j < count; j++)
			{
				for (auto k = 0; k < count; k++)
				{
					double x = (p1[i] < t1.count&&p1[j] < count1&&p1[k] < count1) ? t1.data[p1[i] + p1[j] * count1 + p1[k] * count12] : 0;
					double y = (p2[i] < t2.count&&p2[j] < count2&&p2[k] < count2) ? t2.data[p2[i] + p2[j] * count2 + p2[k] * count22] : 0;
					s += x*y;
				}
			}
		}
		return s;
	}

	// Aлгоритм поиска восхождением к вершине:
	// применение парных замен в перестановках для максимизации корелляции по частотной статистике символов
	template <class T>
	void breakcipher(vector<int> &p1, vector<int> &p2, const T &t1, const T &t2)
	{
		auto current = fitness(p1, p2, t1, t2);

		for (auto loop = true; loop;)
		{
			loop = false;
			for (auto i = 0; i < p1.size() - 1; i++)
			{
				for (auto j = i + 1; j < p1.size(); j++)
				{
					auto t = p1[i];
					p1[i] = p1[j];
					p1[j] = t;
					auto next = fitness(p1, p2, t1, t2);
					if (next > current)
					{
						current = next;
						loop = true;
						clog << '.';
					}
					else
					{
						t = p1[i];
						p1[i] = p1[j];
						p1[j] = t;
					}
				}
			}
			for (auto i = 0; i < p2.size() - 1; i++)
			{
				for (auto j = i + 1; j < p2.size(); j++)
				{
					auto t = p2[i];
					p2[i] = p2[j];
					p2[j] = t;
					auto next = fitness(p1, p2, t1, t2);
					if (next > current)
					{
						current = next;
						loop = true;
						clog << '.';
					}
					else
					{
						t = p2[i];
						p2[i] = p2[j];
						p2[j] = t;
					}
				}
			}
		}
		clog << endl;
	}

	// Применение алгоритма простой замены для файла
	// Символы из стороки key заменяются на символы из строки value
	void replace(const string &sourceFileName, const string &destFileName, const string &key, const string &value)
	{
		ifstream src(sourceFileName);
		ofstream dest(destFileName);
		for (auto ch = src.get(); ch != EOF; ch = src.get())
		{
			auto index = key.find(ch);
			dest << (index < value.length() ? value[index] : ' ');
		}
		src.close();
		dest.close();
	}
}

int main(int argc, char* argv[])
{
	string localeName = "Russian"; // Имя локали
	string cipherFileName = "cipher.txt"; // Имя файла с зашифрованным текстом
	string plainFileName = "plain.txt"; // Имя файла с расшифрованным текстом
	string sampleFileName = "sample.txt"; // Имя файла с образцом открытого текста
	auto useTrigram = false;

	// Обработка параметров командной строки
	for (auto i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "--locale") == 0) localeName = argv[++i]; // Имя локали
		else if (strcmp(argv[i], "--sample") == 0) sampleFileName = argv[++i]; // Имя файла с образцом открытого текста
		else if (strcmp(argv[i], "--cipher") == 0) cipherFileName = argv[++i]; // Имя файла с зашифрованным текстом
		else if (strcmp(argv[i], "--plain") == 0) plainFileName = argv[++i]; // Имя файла с зашифрованным текстом
		else if (strcmp(argv[i], "--3") == 0) useTrigram = true;
	}

	setlocale(LC_ALL, localeName.c_str());

	// Получение частотных статистик триграмм из файлов
	substitution::trigram cipher(cipherFileName);
	substitution::trigram sample(sampleFileName);

	// Получение частотных статистик биграмм из частотных статистик триграмм
	substitution::bigram bicipher(cipher);
	substitution::bigram bisample(sample);

	// Получение частотных статистик символов из частотных статистик биграмм
	substitution::unogram unocipher(bicipher);
	substitution::unogram unosample(bisample);

	// Определение максимального количества символов в алфавитах
	auto count = max(cipher.count, sample.count);

	// Аллокирование и инициализация подстановок
	vector<int> p1(count);
	vector<int> p2(count);

	for (auto i = 0; i < p1.size(); i++) p1[i] = i;
	for (auto i = 0; i < p2.size(); i++) p2[i] = i;

	substitution::breakcipher<substitution::unogram>(p1, p2, unocipher, unosample);
	substitution::breakcipher<substitution::bigram>(p1, p2, bicipher, bisample);
	if (useTrigram) substitution::breakcipher<substitution::trigram>(p1, p2, cipher, sample);

	string key;
	string value;

	for (auto i = 0; i < p1.size(); i++) key.append(1, p1[i] < cipher.count ? cipher.alphabet[p1[i]] : ' ');
	for (auto i = 0; i < p2.size(); i++) value.append(1, p2[i] < sample.count ? sample.alphabet[p2[i]] : ' ');

	substitution::replace(cipherFileName, plainFileName, key, value);

	return 0;
}

