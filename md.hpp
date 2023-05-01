#ifndef MD_HPP_
#define MD_HPP_

#include "variables.hpp"
#include "observer.hpp"

struct MD {
public:
    MD();
    ~MD();
    void run();
private:
    Variables vars;
    Observer obs;
    void makeconf();
    void update_position();
    void periodic();
    void calculate_force();
    void calculate();
};

#endif