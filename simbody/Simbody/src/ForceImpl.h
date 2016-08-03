#ifndef SimTK_SIMBODY_FORCE_IMPL_H_
#define SimTK_SIMBODY_FORCE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Force.h"

namespace SimTK {

// This is what a Force handle points to.
class ForceImpl : public PIMPLImplementation<Force, ForceImpl> {
public:
    ForceImpl() : forces(0), defaultDisabled(false) {}
    ForceImpl(const ForceImpl& clone) {*this = clone;}

    void setDisabledByDefault(bool shouldBeDisabled) 
    {   invalidateTopologyCache();
        defaultDisabled = shouldBeDisabled; }

    bool isDisabledByDefault() const 
    {   return defaultDisabled; }

    virtual ~ForceImpl() {}
    virtual ForceImpl* clone() const = 0;
    virtual bool dependsOnlyOnPositions() const {
        return false;
    }
    ForceIndex getForceIndex() const {return index;}
    const GeneralForceSubsystem& getForceSubsystem() const 
    {   assert(forces); return *forces; }
    void setForceSubsystem(GeneralForceSubsystem& frcsub, ForceIndex ix) {
        forces = &frcsub;
        index  = ix;
    }
    void invalidateTopologyCache() const {
        if (forces) forces->invalidateSubsystemTopologyCache();
    }

    // Every force element must provide the next two methods. Note that 
    // calcForce() must *add in* (+=) its forces to the given arrays.
    virtual void calcForce(const State&         state, 
                           Vector_<SpatialVec>& bodyForces, 
                           Vector_<Vec3>&       particleForces, 
                           Vector&              mobilityForces) const = 0;
    virtual Real calcPotentialEnergy(const State& state) const = 0;

    virtual void realizeTopology    (State& state) const {}
    virtual void realizeModel       (State& state) const {}
    virtual void realizeInstance    (const State& state) const {}
    virtual void realizeTime        (const State& state) const {}
    virtual void realizePosition    (const State& state) const {}
    virtual void realizeVelocity    (const State& state) const {}
    virtual void realizeDynamics    (const State& state) const {}
    virtual void realizeAcceleration(const State& state) const {}
    virtual void realizeReport      (const State& state) const {}
private:
        // CONSTRUCTION
    GeneralForceSubsystem* forces;  // just a reference; no delete on destruction
    ForceIndex             index;

        // TOPOLOGY "STATE"
    // Changing anything here invalidates the topology of the containing
    // force Subsystem and thus of the whole System. 

    // This says whether the Instance-stage "disabled" flag for this force 
    // element should be initially on or off. Most force elements are enabled 
    // by default.
    bool                   defaultDisabled;

        // TOPOLOGY "CACHE"
    // Nothing in the base Impl class.
};



//------------------------------------------------------------------------------
//                    TWO POINT LINEAR SPRING IMPL
//------------------------------------------------------------------------------
class Force::TwoPointLinearSpringImpl : public ForceImpl {
public:
    TwoPointLinearSpringImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real k, Real x0);
    TwoPointLinearSpringImpl* clone() const {
        return new TwoPointLinearSpringImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real k, x0;
};



//------------------------------------------------------------------------------
//                    TWO POINT LINEAR DAMPER IMPL
//------------------------------------------------------------------------------
class Force::TwoPointLinearDamperImpl : public ForceImpl {
public:
    TwoPointLinearDamperImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real damping);
    TwoPointLinearDamperImpl* clone() const {
        return new TwoPointLinearDamperImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real damping;
};



//------------------------------------------------------------------------------
//                    TWO POINT CONSTANT FORCE IMPL
//------------------------------------------------------------------------------
class Force::TwoPointConstantForceImpl : public ForceImpl {
public:
    TwoPointConstantForceImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real force);
    TwoPointConstantForceImpl* clone() const {
        return new TwoPointConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real force;
};



//------------------------------------------------------------------------------
//                    MOBILITY LINEAR SPRING IMPL
//------------------------------------------------------------------------------
class Force::MobilityLinearSpringImpl : public ForceImpl {
public:
    MobilityLinearSpringImpl(const MobilizedBody& body, int coordinate, Real k, Real x0);
    MobilityLinearSpringImpl* clone() const {
        return new MobilityLinearSpringImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    int coordinate;
    Real k, x0;
};



//------------------------------------------------------------------------------
//                       MOBILITY LINEAR DAMPER IMPL
//------------------------------------------------------------------------------
class Force::MobilityLinearDamperImpl : public ForceImpl {
public:
    MobilityLinearDamperImpl(const MobilizedBody& body, int coordinate, Real damping);
    MobilityLinearDamperImpl* clone() const {
        return new MobilityLinearDamperImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    int coordinate;
    Real damping;
};



//------------------------------------------------------------------------------
//                       MOBILITY CONSTANT FORCE IMPL
//------------------------------------------------------------------------------
class Force::MobilityConstantForceImpl : public ForceImpl {
public:
    MobilityConstantForceImpl(const MobilizedBody& body, int coordinate, Real force);
    MobilityConstantForceImpl* clone() const {
        return new MobilityConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    int coordinate;
    Real force;
};



//------------------------------------------------------------------------------
//                           CONSTANT FORCE IMPL
//------------------------------------------------------------------------------
class Force::ConstantForceImpl : public ForceImpl {
public:
    ConstantForceImpl(const MobilizedBody& body, const Vec3& station, const Vec3& force);
    ConstantForceImpl* clone() const {
        return new ConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    Vec3 station, force;
};



//------------------------------------------------------------------------------
//                          CONSTANT TORQUE IMPL
//------------------------------------------------------------------------------
class Force::ConstantTorqueImpl : public ForceImpl {
public:
    ConstantTorqueImpl(const MobilizedBody& body, const Vec3& torque);
    ConstantTorqueImpl* clone() const {
        return new ConstantTorqueImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    Vec3 torque;
};



//------------------------------------------------------------------------------
//                          GLOBAL DAMPER IMPL
//------------------------------------------------------------------------------
class Force::GlobalDamperImpl : public ForceImpl {
public:
    GlobalDamperImpl(const SimbodyMatterSubsystem& matter, Real damping);
    GlobalDamperImpl* clone() const {
        return new GlobalDamperImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    Real damping;
};



//------------------------------------------------------------------------------
//                          UNIFORM GRAVITY IMPL
//------------------------------------------------------------------------------
class Force::UniformGravityImpl : public ForceImpl {
public:
    UniformGravityImpl(const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight);
    UniformGravityImpl* clone() const {
        return new UniformGravityImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    Vec3 getGravity() const {
        return g;
    }
    void setGravity(const Vec3& gravity) {
        g = gravity;
        invalidateTopologyCache();
    }
    Real getZeroHeight() const {
        return zeroHeight;
    }
    void setZeroHeight(Real height) {
        zeroHeight = height;
        invalidateTopologyCache();
    }
private:
    const SimbodyMatterSubsystem& matter;
    Vec3 g;
    Real zeroHeight;
};



//------------------------------------------------------------------------------
//                              CUSTOM IMPL
//------------------------------------------------------------------------------
class Force::CustomImpl : public ForceImpl {
public:
    CustomImpl(Force::Custom::Implementation* implementation);
    CustomImpl* clone() const {
        return new CustomImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return implementation->dependsOnlyOnPositions();
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    ~CustomImpl() {
        delete implementation;
    }
    const Force::Custom::Implementation& getImplementation() const {
        return *implementation;
    }
    Force::Custom::Implementation& updImplementation() {
        return *implementation;
    }
protected:
    void realizeTopology(State& state) const {
        implementation->realizeTopology(state);
    }
    void realizeModel(State& state) const {
        implementation->realizeModel(state);
    }
    void realizeInstance(const State& state) const {
        implementation->realizeInstance(state);
    }
    void realizeTime(const State& state) const {
        implementation->realizeTime(state);
    }
    void realizePosition(const State& state) const {
        implementation->realizePosition(state);
    }
    void realizeVelocity(const State& state) const {
        implementation->realizeVelocity(state);
    }
    void realizeDynamics(const State& state) const {
        implementation->realizeDynamics(state);
    }
    void realizeAcceleration(const State& state) const {
        implementation->realizeAcceleration(state);
    }
    void realizeReport(const State& state) const {
        implementation->realizeReport(state);
    }
private:
    Force::Custom::Implementation* implementation;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_IMPL_H_
