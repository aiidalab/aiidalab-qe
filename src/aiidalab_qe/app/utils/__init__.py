from importlib.metadata import distributions


def print_error(entry_point, e):
    print(f"\033[91mFailed to load plugin entry point {entry_point.name}.\033[0m")
    print(
        "\033[93mThis may be due to compatibility issues with the current QEApp version.\033[0m"
    )
    print("\033[93mPlease contact the plugin author for further assistance.\033[0m")
    print(
        "\033[93mThus, the plugin will not be available. However, you can still use the rest of the app.\033[0m"
    )
    print(f"\033[91mError message: {e}\033[0m\n")


# load entry points
def get_entries(entry_point_name="aiidalab_qe.properties"):
    from importlib.metadata import entry_points

    entries = {}
    for entry_point in entry_points().get(entry_point_name, []):
        try:
            # Attempt to load the entry point
            loaded_entry_point = entry_point.load()
            entries[entry_point.name] = loaded_entry_point
        except Exception as e:
            print_error(entry_point, e)

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
