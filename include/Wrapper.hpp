#pragma once
#include <memory>
#include <vector>
#include <utility>

// Wrapped by Jinkela and Eric 2023.11.28

class Pin;
class Net;
class Module;
class Placement;

namespace wrapper
{
    class Pin
    {
    public:
        Pin(::Pin &_origin);
        ~Pin();

        double x();
        double y();
        unsigned moduleId();
        unsigned netId();

        ::Pin &getOrigin();
        const ::Pin &getOrigin() const;

    private:
        ::Pin &_origin;
    };

    class Net
    {
    public:
        Net(::Net &_origin);
        ~Net();

        unsigned numPins();
        Pin pin(unsigned index); // index: 0 ~ numPins()-1

        ::Net &getOrigin();
        const ::Net &getOrigin() const;

    private:
        ::Net &_origin;
    };

    class Module
    {
    public:
        Module(::Module &_origin);
        ~Module();

        const char *name();
        double x();
        double y();
        double width();
        double height();
        bool isFixed();
        const char *orientString();

        double centerX();
        double centerY();
        double area();

        unsigned numPins();
        Pin pin(unsigned index); // index: 0 ~ numPins()-1

        void setPosition(double x, double y);
        void setCenterPosition(double x, double y);

        ::Module &getOrigin();
        const ::Module &getOrigin() const;

    private:
        ::Module &_origin;
    };

    class Placement
    {
    public:
        Placement();
        ~Placement();

        const char *name();
        double boundryTop();
        double boundryLeft();
        double boundryBottom();
        double boundryRight();

        unsigned numModules();
        unsigned numNets();
        unsigned numPins();

        Module module(unsigned moduleId);
        Net net(unsigned netId);
        Pin pin(unsigned pinId);

        double computeHpwl();

        void readBookshelfFormat(const char *const filePathName, const char *const plFileName);
        void outputBookshelfFormat(const char *const filePathName);

        ::Placement &getOrigin();
        const ::Placement &getOrigin() const; 

    private:
        std::unique_ptr<::Placement> _origin;
 
    };
}