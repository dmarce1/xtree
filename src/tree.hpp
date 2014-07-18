/*
 * tree.hpp
 *
 *  Created on: May 24, 2014
 *      Author: dmarce1
 */

#ifndef TREE_HPP_
#define TREE_HPP_

namespace xtree {

template<typename Derived, int Ndim>
class tree: public hpx::components::managed_component_base<tree<Derived, Ndim>, hpx::components::detail::this_type,
		hpx::traits::construct_with_back_ptr> {
public:
	static constexpr int Nbranch = 2;
	static constexpr int Nneighbor = pow_<3, Ndim>::value;
	const char* name = "tree";
	const char* silo_name = "silo_output";
private:
	using base_type = hpx::components::managed_component_base<tree<Derived,Ndim>, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr>;
	using component_type = hpx::components::managed_component<tree<Derived,Ndim>>;
	using silo_output_type = silo_output<Ndim>;
	hpx::id_type silo_gid;
	hpx::id_type this_gid;
	hpx::id_type root_node_gid;
	std::set<Derived*> nodes;
	tree<Derived, Ndim>* this_ptr;
	std::array<hpx::id_type, Nbranch> child_gids;
	load_balancer* load_balancer_ptr;
	hpx::id_type load_balancer_gid;
	mutable hpx::lcos::local::mutex dir_lock;
public:

	tree() {
		assert(false);
	}
	tree(component_type* back_ptr) :
			base_type(back_ptr) {
		static bool initialized = false;
		assert(!initialized);
		initialized = true;
		int my_id, id_cnt;
		auto localities = hpx::find_all_localities();
		std::vector < hpx::future < hpx::id_type >> futures(Nbranch);
//		printf("1\n");
		my_id = hpx::get_locality_id();
		id_cnt = localities.size();
		this_gid = base_type::get_gid();
		this_ptr = this;
		hpx::register_id_with_basename(name, this_gid, my_id).get();
//		printf("2\n");
		for (int i = 0; i < Nbranch; i++) {
			int j = my_id * Nbranch + i + 1;
			if (j < localities.size()) {
				futures[i] = hpx::new_<tree<Derived, Ndim>>(localities[j]);
			} else {
				futures[i] = hpx::make_ready_future(hpx::invalid_id);
			}
		}
//		printf("3\n");
		load_balancer_gid = hpx::new_ < load_balancer > (hpx::find_here()).get();
//		printf("3.4\n");
		load_balancer_ptr = (hpx::async < load_balancer::action_get_ptr > (load_balancer_gid)).get();
	//	printf("3.6\n");
	//	printf("3.7\n");
		if (my_id == 0) {
			hpx::future < hpx::id_type > fut1;
	//		printf("silonewin\n");
			fut1 = hpx::new_ < silo_output_type > (hpx::find_here());
	//		printf("silonewout\n");
			silo_gid = fut1.get();
	//		printf("silonewout2\n");
			hpx::register_id_with_basename(silo_name, silo_gid, 0).get();
		} else {
			silo_gid = (hpx::find_id_from_basename(silo_name, 0)).get();
		}
		for (int i = 0; i < Nbranch; i++) {
			child_gids[i] = futures[i].get();
		}
//		printf("4\n");
		root_node_gid = hpx::invalid_id;
	}

	void place_root() {
		hpx::future < hpx::id_type > fut2;
		std::vector < hpx::id_type > neighbors( Nneighbor);
		std::fill(neighbors.begin(), neighbors.end(), hpx::invalid_id);
		fut2 = new_node(location<Ndim>(), hpx::invalid_id, neighbors, 0);
		root_node_gid = fut2.get();
	}

	hpx::future<Derived*> get_root() {
		if (root_node_gid != hpx::invalid_id) {
			return hpx::async<typename node<Derived, Ndim>::action_get_this>(root_node_gid);
		} else {
			return hpx::make_ready_future<Derived*>(nullptr);
		}
	}

	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		assert(false);
	}

	virtual ~tree() {
		if (root_node_gid != hpx::invalid_id) {
			hpx::async<typename node<Derived, Ndim>::action_debranch>(root_node_gid).get();
		}
	}

	tree* get_this() {
		return this;
	}

	hpx::id_type get_new_node(const location<Ndim>& _loc, hpx::id_type _parent_id,
			const std::vector<hpx::id_type>& _neighbors, int subcyc) {
		hpx::shared_future < hpx::id_type > id_future;
		auto fut0 = hpx::new_ < Derived > (hpx::find_here());
		id_future = fut0.share();
		auto fut1 =
				id_future.then(
						hpx::util::unwrapped(
								[=](hpx::id_type id) {
									return hpx::async<typename node<Derived,Ndim>::action_initialize>(id, _loc, hpx::util::make_tuple(_parent_id, std::move(_neighbors)), this);
								}));
		return fut1.then(hpx::util::unwrapped([this,id_future](Derived* ptr) {
			dir_lock.lock();
			auto test = nodes.insert(ptr);
			dir_lock.unlock();
			assert(test.second);
			return id_future.get();
		})).get();
	}

	hpx::future<hpx::id_type> new_node(const location<Ndim>& _loc, hpx::id_type _parent_id,
			const std::vector<hpx::id_type>& _neighbors, int subcyc) {
		auto proc_num = load_balancer_ptr->increment_load().get();
		auto gid = hpx::find_id_from_basename(name, proc_num).get();
		return hpx::async < action_get_new_node > (gid, _loc, _parent_id, _neighbors, subcyc);
	}

	void delete_node(Derived* ptr) {
		load_balancer_ptr->decrement_load();
		dir_lock.lock();
		auto iter = nodes.find(ptr);
		assert(iter != nodes.end());
		nodes.erase(iter);
		dir_lock.unlock();
	}
	void output() const {
		for (int i = 0; i < Nbranch; i++) {
			if (child_gids[i] != hpx::invalid_id) {
				hpx::apply < action_output > (child_gids[i]);
			}
		}
		dir_lock.lock();

		std::size_t leaf_cnt = 0;
		for (auto i = nodes.begin(); i != nodes.end(); ++i) {
			if ((*i)->is_terminal()) {
				leaf_cnt++;
			}
		}
		printf("Leaf cnt = %li\n", leaf_cnt);
		std::vector<typename silo_output_type::zone> zones(leaf_cnt * Derived::Size);
		auto j = zones.begin();
		for (auto i = nodes.begin(); i != nodes.end(); ++i) {
			if ((*i)->is_terminal()) {
				const auto these_zones = (*i)->get_output_zones();
				for (auto k = these_zones.begin(); k != these_zones.end(); ++k) {
					*j = *k;
					++j;
				}
			}
		}

		auto fut = hpx::async<typename silo_output_type::action_send_zones_to_silo>(silo_gid, hpx::get_locality_id(),
				zones);
		fut.get();
		dir_lock.unlock();

	}

	HPX_DEFINE_COMPONENT_ACTION_TPL( tree,get_new_node,action_get_new_node ); //
	HPX_DEFINE_COMPONENT_ACTION_TPL( tree,get_this,action_get_this ); //
	HPX_DEFINE_COMPONENT_ACTION_TPL( tree,place_root, action_place_root); //
	HPX_DEFINE_COMPONENT_ACTION_TPL( tree,output,action_output );
//

};


}

#endif /* TREE_HPP_ */
