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
    from importlib_metadata import entry_points

    entries = {}
    for entry_point in entry_points(group=entry_point_name):
        try:
            # Attempt to load the entry point
            if entry_point.name in entries:
                continue
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
        raise ValueError(f"Package '{package_name}' not found.")
    # Retrieve all entry points associated with this distribution
    if dist.entry_points:
        for ep in dist.entry_points:
            if ep.group == group:
                entry_points_list.append(ep)
    else:
        print(f"No entry points found for package '{package_name}'.")
    return entry_points_list


def test_plugin_functionality(plugin_name):
    """Test the functionality of the plugin.
    1) loading all entry points.
    2) check if the plugin use correct QEAppComputationalResourcesWidget
    """
    from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget

    try:
        eps = get_entry_points_for_package(plugin_name)
        # check if we can load all entry points
        for ep in eps:
            loaded_ep = ep.load()
            # check if the code uses the correct widget
            for name, code in loaded_ep.get("code", {}).items():
                if not isinstance(code, QEAppComputationalResourcesWidget):
                    return (
                        False,
                        f"\nPlugin {plugin_name} code {name} must use QEAppComputationalResourcesWidget class",
                    )
    except Exception as e:
        return False, f"Failed to get entry points for package {plugin_name}: {e}"
    return True, ""
