#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <boost/random.hpp>

class RandomGenerator{
public:
    RandomGenerator(int seed, int len);
    int unidir();
    int unidir_loop();
    int unimove();
    int unimove2();
    int unisitering(); //now can be easily modified for a better choice of random sites (to have less uncounted moves)
    double disReal();

private:
    int len; //polymer length
    boost::random::mt19937_64 gen;
    boost::random::uniform_01<double> unireal;
    boost::uniform_int<int> uniint01{0,1};
    boost::uniform_int<int> uniint02{0,2};
    boost::uniform_int<int> uniint15{1,5};
};

RandomGenerator::RandomGenerator(int seed, int len) {
    this -> len = len;
    gen.seed(seed);
}
int RandomGenerator::unidir() {
    return uniint02(gen);
}
int RandomGenerator::unidir_loop() {
    return uniint15(gen);
}
int RandomGenerator::unimove() {
    return uniint02(gen);
}
int RandomGenerator::unimove2() {
    return uniint01(gen);
}
int RandomGenerator::unisitering() {
    return std::uniform_int_distribution<int>{0, len-1}(gen);
}
double RandomGenerator::disReal() {
    return unireal(gen);
}

#endif 
