from importlib.metadata import distributions


# load entry points
def get_entries(entry_point_name="aiidalab_qe.properties"):
    from importlib.metadata import entry_points

    entries = {}
    for entry_point in entry_points().get(entry_point_name, []):
        entries[entry_point.name] = entry_point.load()

    return entries


def get_entry_items(entry_point_name, item_name="outline"):
    entries = get_entries(entry_point_name)
    return {
        name: entry_point.get(item_name)
        for name, entry_point in entries.items()
        if entry_point.get(item_name, False)
    }


def get_entry_points_for_package(
    package_name: str, group: str = "aiidalab_qe.properties"
):
    """Find the entry points for the specified package"""
    entry_points_list = []

    dist = next(
        (d for d in distributions() if d.metadata["Name"] == package_name), None
    )
    if not dist:
        print(f"No distribution found for package '{package_name}'.")
        return entry_points_list

    # Retrieve all entry points associated with this distribution
    if dist.entry_points:
        for ep in dist.entry_points:
            if ep.group == group:
                entry_points_list.append(ep)
    else:
        print(f"No entry points found for package '{package_name}'.")
    return entry_points_list


def test_plugin_functionality(plugin_name):
    """Test the functionality of the plugin by loading all entry points."""
    try:
        eps = get_entry_points_for_package(plugin_name)
        # check if we can load all entry points
        try:
            for ep in eps:
                ep.load()
        except Exception as e:
            return False, f"Failed to load entry point {ep.name}: {e}"
    except Exception as e:
        return False, f"Failed to get entry points for package {plugin_name}: {e}"
    return True, ""
