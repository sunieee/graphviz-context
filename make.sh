make -j8
sudo make install
cat tests/test/test.dot | dot -v -Tsvg > tests/test/test.svg 2> tests/test/test.log
cat tests/test/test1.dot | dot -v -Tsvg > tests/test/test1.svg 2> tests/test/test1.log