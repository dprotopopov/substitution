﻿# Реализация атаки на шифр простой замены
Методы шифрования заменой (подстановкой) основаны на том, что символы исходного текста, обычно разделенные на блоки и записанные в одном алфавите, заменяются одним или несколькими символами другого алфавита в соответствии с принятым правилом преобразования.
В простом случае при шифровании заменой символы исходного текста, записанные в одном алфавите A, заменяются символами другого алфавита B в соответствии с принятым правилом преобразования.

## Одноалфавитная замена
Одним из важных подклассов методов замены являются одноалфавитные (или моноалфавитные) подстановки, в которых устанавливается однозначное соответствие между каждым знаком ai исходного алфавита сообщений A и соответствующим знаком ei зашифрованного текста E. Одноалфавитная подстановка иногда называется также простой заменой, так как является самым простым шифром замены.
Примером одноалфавитной замены является шифр Цезаря.

# Криптоанализ шифра простой замены
Пространство ключей шифра простой замены огромно и равно количеству перестановок используемого алфавита. Так для английского языка это число составляет 26! = 2^88. Разумеется наивный перебор всех возможных ключей дело безнадежное и для взлома потребуется более утонченная техника, такая как поиск восхождением к вершине:
1. Выбирается случайная перестановка букв — основной ключ. Шифртекст расшифровывается с помощью основного ключа. Для получившегося текста вычисляется коэффициент, характеризующий вероятность принадлежности к естественному языку.
2. Основной ключ подвергается небольшим изменениям (перестановка двух произвольно выбранных букв). Производится расшифровка и вычисляется коэффициент полученного текста.
3. Если коэффициент выше сохраненного значения, то основной ключ заменяется на модифицированный вариант.
4. Шаги 2-3 повторяются пока коэффициент не станет постоянным.

# Алгоритм вскрытия шифра простой замены
Применяется алгоритм поиска восхождением к вершине для вычисления максимума корелляции между частотными статистиками образца открытого текста и шифротекста.
Символы алфавитов A (открытого текста) и B (шифротекста) могут быть определены из образца открытого текста и шифротекста.
Функция дешифрования D:B->A может быть представлена в виде композиции двух подстановок D=(P1^-1)*P2, где P1:N->B, P2:N->A.
Подстановки P1 и P2 являются точкой максимума функции корелляции, подсчитанной по частотным статистикам образца открытого текста и шифротекста.
Для вычисления подстановок P1 и P2 будем последовательно заменять у них пары аргументов если такая замена увеличивает корелляцию.
Будем повторять это действие до тех пор, пока происходит увеличение корелляции.
Полученные таким образом подстановки P1 и P2 являются точкой максимума функции корелляции, а значит может быть вычисленна функция дешифрования D:B->A.

# Ссылки

1. Криптоанализ классических шифров на сайте http://practicalcryptography.com
2. Описание поиска восхождением к вершине на https://en.wikipedia.org/wiki/Hill_climbing