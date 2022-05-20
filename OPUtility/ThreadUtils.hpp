#pragma once

#include <mutex>
#include <thread>
#include <map>
#include <ranges>

/**
 * \brief Synchronized counter implementation, starting from n.
 */
class SynchronizedCounter
{
	int c_;
	mutable std::mutex mc_;

public:
	explicit SynchronizedCounter(const int& n = 0) : c_(n) {}

	SynchronizedCounter& operator++();

	SynchronizedCounter& increment();

	int incrementAndGet();
};

inline SynchronizedCounter& SynchronizedCounter::operator++()
{
	std::lock_guard g(mc_);
	this->c_ += 1;
	return *this;
}

inline SynchronizedCounter& SynchronizedCounter::increment()
{
	return ++(*this);
}

inline int SynchronizedCounter::incrementAndGet()
{
	++(*this);
	return this->c_;
}


/**
 * \brief Wrapper JThread with name.
 */
class NamedJThread : public std::jthread
{
	inline static SynchronizedCounter counter_;

	std::string name_;

public:
	template <typename Fn, typename... Args>
	explicit NamedJThread(Fn&& func, Args&&... args) : jthread(func, args...), name_("Executor " + std::to_string(counter_.incrementAndGet())) {}

	template <typename Fn, typename... Args>
	explicit NamedJThread(Fn&& func, Args&&... args, const std::string& name) : jthread(func, args...), name_(name) {}

	[[nodiscard]] std::string getName() const;
	void setName(const std::string& name);
};

inline std::string NamedJThread::getName() const
{
	return this->name_;
}

inline void NamedJThread::setName(const std::string& name)
{
	this->name_ = name;
}


/**
 * \brief Arbitary thread container of elements derived from std::jthread.
 */
template <typename ThreadType>
class ExecutorContainer final
{
	static_assert(std::is_base_of_v<std::jthread, ThreadType>, "ThreadType must inherit from std::jthread.");

	std::map<std::thread::id, std::shared_ptr<ThreadType>> executors_;
	mutable std::mutex m1Executors_;
	mutable std::mutex m2Executors_;

public:
	template <typename Fn, typename... Args>
	void createExecutor(Fn&& func, Args&&... args);

	[[nodiscard]] std::optional<const ThreadType*> getById(const std::thread::id& id) const;

	[[nodiscard]] std::vector<const ThreadType*> getExecutors() const;

	void waitForExecutors() const;
	void clearExecutors();
};

template <typename ThreadType>
template <typename Fn, typename... Args>
void ExecutorContainer<ThreadType>::createExecutor(Fn&& func, Args&&... args)
{
	std::lock_guard g1(m1Executors_);
	std::lock_guard g2(m2Executors_);
	auto th(std::make_shared<ThreadType>(func, args...));
	std::thread::id currId(th->get_id());
	executors_.insert(std::make_pair(currId, th));
}

template <typename ThreadType>
std::optional<const ThreadType*> ExecutorContainer<ThreadType>::getById(const std::thread::id& id) const
{
	static_assert(std::is_base_of_v<std::jthread, ThreadType>, "ThreadType must inherit from std::jthread.");

	std::lock_guard g1(m1Executors_);
	if (executors_.contains(id)) return std::optional(executors_.at(id).get());
	return std::nullopt;
}

template <typename ThreadType>
std::vector<const ThreadType*> ExecutorContainer<ThreadType>::getExecutors() const
{
	static_assert(std::is_base_of_v<std::jthread, ThreadType>, "ThreadType must inherit from std::jthread.");

	std::lock_guard g1(m1Executors_);
	std::vector<const ThreadType*> lexecutors;

	for (const auto& executor : executors_ | std::views::values)
	{
		lexecutors.push_back(executor.get());
	}

	return lexecutors;
}

template <typename ThreadType>
void ExecutorContainer<ThreadType>::waitForExecutors() const
{
	static_assert(std::is_base_of_v<std::jthread, ThreadType>, "ThreadType must inherit from std::jthread.");

	std::lock_guard g2(m2Executors_);
	for (const auto& executor : executors_ | std::views::values)
	{
		executor->join();
	}
}

template <typename ThreadType>
void ExecutorContainer<ThreadType>::clearExecutors()
{
	static_assert(std::is_base_of_v<std::jthread, ThreadType>, "ThreadType must inherit from std::jthread.");

	std::lock_guard g1(m1Executors_);
	std::lock_guard g2(m2Executors_);

	executors_.clear();
}
