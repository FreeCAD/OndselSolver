// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace MbD {
// A run-local event scheduler. Geometry/measurement and action adapters belong
// to the caller. Probing dense output never commits an action or changes state.
class DynamicEvents {
public:
    struct Event {
        enum Trigger { Time, Threshold, After } trigger = Time;
        std::string name;
        double time = 0, threshold = 0, hysteresis = 0, delay = 0;
        size_t predecessor = 0;
        bool rising = true, repeat = false, fireInitially = false;
        std::function<double()> measure;
        std::function<void(double)> action;
        bool fired = false, armed = true;
        double lastFire = -std::numeric_limits<double>::infinity();
        std::deque<double> scheduledTimes;
    };
    struct Firing { size_t event; double time; };
    std::vector<Event> events;
    std::vector<Firing> firings;
    std::function<void()> prepare;
    double tolerance = 1e-9;

    bool initialize(double time) {
        pending.clear();
        firings.clear();
        beforeActions.clear();
        if (!std::isfinite(tolerance) || tolerance <= 0)
            throw std::invalid_argument("Event tolerance must be finite and positive");
        for (auto& event : events) {
            if (!event.action || !std::isfinite(event.time) || !std::isfinite(event.delay)
                || event.delay < 0 || !std::isfinite(event.threshold)
                || !std::isfinite(event.hysteresis) || event.hysteresis < 0
                || (event.trigger == Event::Threshold && !event.measure)
                || (event.trigger == Event::After && event.predecessor >= events.size()))
                throw std::invalid_argument("Invalid event: " + event.name);
            event.fired = false;
            event.armed = true;
            event.lastFire = -std::numeric_limits<double>::infinity();
            event.scheduledTimes.clear();
        }
        for (size_t i = 0; i < events.size(); ++i) {
            auto& event = events[i];
            if (event.trigger == Event::Threshold) {
                const double v = measure(event);
                event.armed = !satisfied(event, v);
                if (event.fireInitially && satisfied(event, v)) pending.push_back(i);
            }
            else if (event.trigger == Event::Time && std::abs(event.time-time) <= tolerance) {
                pending.push_back(i);
            }
        }
        if (pending.empty()) return false;
        apply(time);
        return true;
    }

    // The step is converged, but not yet committed. At most eight equal probe
    // intervals are searched; fast oscillatory signals require a smaller hmax.
    double locate(double previous, double current, const std::function<void(double)>& sample) {
        pending.clear();
        double earliest = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < events.size(); ++i) {
            const auto& event = events[i];
            if (event.fired && !event.repeat) continue;
            double candidate = std::numeric_limits<double>::infinity();
            if (event.trigger != Event::Threshold) {
                candidate = scheduled(event);
                if (candidate <= previous+tolerance || candidate > current+tolerance) continue;
            }
            else {
                bool armed = event.armed;
                sample(previous);
                double left = previous, vleft = measure(event);
                for (int probe = 1; probe <= 8; ++probe) {
                    const double right = previous+(current-previous)*probe/8;
                    sample(right);
                    const double vright = measure(event);
                    if (!armed && rearmed(event, vleft)) armed = true;
                    if (armed && !satisfied(event, vleft) && satisfied(event, vright)) {
                        double lo = left, hi = right;
                        for (int iter = 0; iter < 60 && hi-lo > tolerance; ++iter) {
                            const double mid = (lo+hi)/2;
                            sample(mid);
                            if (satisfied(event, measure(event))) hi = mid;
                            else lo = mid;
                        }
                        candidate = hi;
                        break;
                    }
                    left = right;
                    vleft = vright;
                }
            }
            if (!std::isfinite(candidate)) continue;
            if (candidate < earliest-tolerance) {
                earliest = candidate;
                pending.clear();
            }
            if (std::abs(candidate-earliest) <= tolerance) pending.push_back(i);
        }
        sample(current);
        return earliest;
    }

    void commit() {
        for (auto& event : events) {
            if (event.trigger == Event::Threshold && !event.armed && rearmed(event, measure(event)))
                event.armed = true;
        }
    }

    // Rearming can occur inside a step, even if its endpoint has already
    // returned into the hysteresis band. Commit only the accepted interval.
    void commit(double previous, double current, const std::function<void(double)>& sample) {
        for (int probe = 0; probe <= 8; ++probe) {
            sample(previous + (current-previous)*probe/8);
            commit();
        }
    }

    void apply(double time) {
        beforeActions.clear();
        for (const auto& event : events)
            beforeActions.push_back(event.trigger == Event::Threshold ? measure(event) : 0.0);
        auto batch = pending;
        pending.clear();
        for (size_t pass = 0; !batch.empty(); ++pass) {
            if (pass > events.size() || firings.size()+batch.size() > 10000)
                throw std::runtime_error("Event cascade exceeds the safe execution limit");
            for (auto index : batch) {
                auto& event = events.at(index);
                event.action(time);
                event.fired = true;
                event.armed = false;
                event.lastFire = time;
                if (event.trigger == Event::After && !event.scheduledTimes.empty())
                    event.scheduledTimes.pop_front();
                firings.push_back({index, time});
                // Preserve each pending occurrence. A later predecessor firing
                // must not reset an already-running delay.
                for (auto& dependent : events) {
                    if (dependent.trigger == Event::After && dependent.predecessor == index
                        && (dependent.repeat || (!dependent.fired && dependent.scheduledTimes.empty())))
                        dependent.scheduledTimes.push_back(time + dependent.delay);
                }
            }
            batch.clear();
            for (size_t i = 0; i < events.size(); ++i) {
                const auto& event = events[i];
                if (event.trigger == Event::After && (!event.fired || event.repeat)
                    && std::abs(scheduled(event)-time) <= tolerance)
                    batch.push_back(i);
            }
        }
        commit();
    }

    // Constraint activation may change velocity instantaneously during IC.
    // Detect those crossings after projection, before advancing physical time.
    bool settle(double time) {
        pending.clear();
        if (beforeActions.size() != events.size()) return false;
        for (size_t i = 0; i < events.size(); ++i) {
            const auto& event = events[i];
            if (event.trigger == Event::Threshold && event.armed
                && (!event.fired || event.repeat) && !satisfied(event, beforeActions[i])
                && satisfied(event, measure(event))) pending.push_back(i);
        }
        beforeActions.clear();
        if (pending.empty()) { commit(); return false; }
        apply(time);
        return true;
    }

private:
    std::vector<size_t> pending;
    std::vector<double> beforeActions;
    static double measure(const Event& event) {
        const double result = event.measure();
        if (!std::isfinite(result)) throw std::runtime_error("Non-finite event measurement: "+event.name);
        return result;
    }
    static bool satisfied(const Event& e, double value) {
        return e.rising ? value >= e.threshold : value <= e.threshold;
    }
    static bool rearmed(const Event& e, double value) {
        return e.rising ? value < e.threshold-e.hysteresis : value > e.threshold+e.hysteresis;
    }
    double scheduled(const Event& event) const {
        if (event.trigger == Event::Time) return event.fired ? std::numeric_limits<double>::infinity() : event.time;
        return event.scheduledTimes.empty() ? std::numeric_limits<double>::infinity()
                                           : event.scheduledTimes.front();
    }
};
}
