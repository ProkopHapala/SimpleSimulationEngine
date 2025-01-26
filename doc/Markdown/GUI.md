# GUI System Documentation

This document provides a comprehensive guide to the GUI system implemented in `GUI.h` and `GUI.cpp`, with usage examples from `MolGUI.h`. It serves as both a reference and a tutorial for using this library in other projects.

## Design Philosophy

The GUI system is designed to be a flexible and extensible framework for creating interactive user interfaces within OpenGL applications. It uses a panel-based approach, where different UI elements are represented as panels that can be arranged and interacted with. The system is event-driven, responding to user input such as mouse clicks and key presses.

## Key Classes

### `GUI`

The `GUI` class is the main container for all GUI elements. It manages a list of `GUIAbstractPanel` objects and handles global events.

-   `addPanel(GUIAbstractPanel* panel)`: Adds a panel to the GUI.
-   `onEvent(int mouseX, int mouseY, const SDL_Event& event)`: Handles SDL events and dispatches them to the appropriate panels.
-   `draw()`: Renders all visible panels.
-   `layoutRow(int xmin, int ymin, int xspace=0)`: Arranges panels in a row.
-   `clear(int i0=0)`: Clears all panels from the GUI.

### `GUIAbstractPanel`

This is the base class for all GUI panels. It provides common functionality for rendering, event handling, and layout.

-   `xmin`, `xmax`, `ymin`, `ymax`: Defines the panel's bounding box.
-   `visible`, `disabled`, `redraw`: Flags for controlling panel state.
-   `bgColor`, `textColor`: Colors for the panel background and text.
-   `caption`: The panel's title.
-   `command`: A callback function to be executed when the panel is activated.
-   `initPanel(const std::string& caption_, int xmin_, int ymin_, int xmax_, int ymax_)`: Initializes the panel's properties.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.
-   `onKeyDown(const SDL_Event& e, GUI& gui)`: Handles key press events within the panel.
-   `onText(const SDL_Event& e, GUI& gui)`: Handles text input events within the panel.
-   `view()`: Renders the panel.
-   `render()`: Renders the panel's content.
-   `draw()`: Calls `tryRender()` and `view()` to render the panel.
-   `check(int x, int y)`: Checks if a point is within the panel's bounds.
-   `toRelative(int& x, int& y)`: Converts screen coordinates to panel-relative coordinates.
-   `setCommand(const std::function<void(GUIAbstractPanel* panel)>& command_)`: Sets the panel's command callback.

### `GUIPanel`

A concrete panel class that implements a slider or button.

-   `isSlider`, `isButton`, `bCmdOnSlider`, `isInt`, `viewVal`, `valBelow`: Flags for controlling panel behavior.
-   `barColor`: Color of the slider bar.
-   `value`: The current value of the slider.
-   `vmin`, `vmax`: The minimum and maximum values of the slider.
-   `master`: A pointer to a variable that the slider controls.
-   `view()`: Renders the slider or button.
-   `render()`: Renders the panel's content.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.
-   `onKeyDown(const SDL_Event& e, GUI& gui)`: Handles key press events within the panel.
-   `onText(const SDL_Event& e, GUI& gui)`: Handles text input events within the panel.
-   `x2val(float x)`: Converts a horizontal position to a slider value.
-   `val2x(double val)`: Converts a slider value to a horizontal position.
-   `syncRead()`: Reads the value from the master variable.
-   `syncWrite()`: Writes the value to the master variable.
-   `setRange(float vmin_, float vmax_)`: Sets the slider's range.
-   `setValue(float val_)`: Sets the slider's value.

### `MultiPanel`

A panel that contains multiple sub-panels, typically `GUIPanel` objects.

-   `nsubs`: The number of sub-panels.
-   `subs`: A vector of sub-panels.
-   `opened`: A flag indicating whether the sub-panels are visible.
-   `dy`: The vertical spacing between sub-panels.
-   `initMulti(const std::string& caption, int xmin_, int ymin_, int xmax_, int dy_, int nsubs_, bool isSlider, bool isButton, bool isInt, bool viewVal, bool bCmdOnSlider)`: Initializes the multi-panel.
-   `addPanel(const std::string& label, Vec3d vals, bool isSlider_=true, bool isButton_=true, bool isInt_=false, bool viewVal_=true, bool bCmdOnSlider_=false)`: Adds a sub-panel.
-   `open()`: Opens the multi-panel, making sub-panels visible.
-   `close()`: Closes the multi-panel, hiding sub-panels.
-   `toggleOpen()`: Toggles the open/close state of the multi-panel.
-   `view()`: Renders the multi-panel and its sub-panels.
-   `render()`: Renders the multi-panel's content.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.

### `CheckBoxList`

A panel that displays a list of checkboxes.

-   `boxes`: A vector of `CheckBox` objects.
-   `dy`: The vertical spacing between checkboxes.
-   `checkColor`: The color of the checkbox.
-   `addBox(std::string label, bool* ptr)`: Adds a checkbox to the list.
-   `initCheckBoxList(int xmin_, int ymin_, int xmax_, int dy=fontSizeDef*2)`: Initializes the checkbox list.
-   `update()`: Updates the checkbox list.
-   `view()`: Renders the checkbox list.
-   `render()`: Renders the checkbox list's content.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.
-   `syncRead()`: Reads the values from the master variables.
-   `syncWrite()`: Writes the values to the master variables.

### `ScisorBox`

A panel that defines a scissor rectangle for OpenGL rendering.

-   `apply()`: Applies the scissor rectangle.
-   `initScisor(const std::string& caption, int xmin_, int ymin_, int xmax_, int ymax_)`: Initializes the scisor box.
-   `render()`: Renders the scisor box's content.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.

### `CommandList`

A panel that displays a list of commands with key bindings.

-   `commands`: A pointer to a `Commander` object.
-   `commandDispatch`: A callback function to be executed when a command is selected.
-   `dy`: The vertical spacing between commands.
-   `modColor`: The color of the selected command.
-   `initCommandList(int xmin_, int ymin_, int xmax_, int dy=fontSizeDef*2)`: Initializes the command list.
-   `getKeyb(int key)`: Handles key press events within the panel.
-   `update()`: Updates the command list.
-   `view()`: Renders the command list.
-   `render()`: Renders the command list's content.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.

### `DropDownList`

A panel that displays a drop-down list of items.

-   `bOpened`: A flag indicating whether the list is open.
-   `nSlots`: The number of visible items in the list.
-   `iSelected`: The index of the selected item.
-   `labels`: A vector of item labels.
-   `addItem(const std::string& label)`: Adds an item to the list.
-   `initList(const std::string& caption, int xmin_, int ymin_, int xmax_, int nSlots_)`: Initializes the drop-down list.
-   `selectedToStr(char* str)`: Gets the string representation of the selected item.
-   `open()`: Opens the drop-down list.
-   `close()`: Closes the drop-down list.
-   `render()`: Renders the drop-down list's content.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.

### `TreeView`

A panel that displays a tree-like structure.

-   `iItem0`: The index of the first visible item.
-   `iSelected`: The index of the selected item.
-   `nSlots`: The number of visible items.
-   `root`: The root node of the tree.
-   `lines`: A vector of visible tree nodes.
-   `updateLines(TreeViewTree& node, int level)`: Updates the list of visible tree nodes.
-   `initTreeView(const std::string& caption, int xmin_, int ymin_, int xmax_, int nSlots_)`: Initializes the tree view.
-   `view()`: Renders the tree view.
-   `render()`: Renders the tree view's content.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.

### `TableView`

A panel that displays a table of data.

-   `table`: A pointer to a `Table` object.
-   `i0`, `j0`: The starting row and column indices.
-   `imax`, `jmax`: The ending row and column indices.
-   `i`, `j`: The current row and column indices.
-   `nchs`: A vector of column widths.
-   `xs`: A vector of column x-positions.
-   `input`: A pointer to a `GUITextInput` object for editing cells.
-   `initTableView(Table* table_, const std::string& caption_, int xmin_, int ymin_, int i0_, int j0_, int imax_, int jmax_)`: Initializes the table view.
-   `render()`: Renders the table view's content.
-   `onKeyDown(const SDL_Event& e, GUI& gui)`: Handles key press events within the panel.
-   `onText(const SDL_Event& e, GUI& gui)`: Handles text input events within the panel.
-   `onMouse(int x, int y, const SDL_Event& event, GUI& gui)`: Handles mouse events within the panel.
-   `view()`: Renders the table view.

### `GUITextInput`

A class that implements a text input field.

-   `curPos`: The current cursor position.
-   `inputText`: The current text in the input field.
-   `isNumber`, `checkRange`: Flags for controlling input behavior.
-   `vmin`, `vmax`: The minimum and maximum values for numeric input.
-   `value`: The current numeric value.
-   `num_op`: The current numeric operation.
-   `modified`, `entered`: Flags for tracking input state.
-   `applyVal(float f)`: Applies a numeric operation to the current value.
-   `view3D(const Vec3d& pos, int fontTex, float textSize)`: Renders the text input in 3D.
-   `viewHUD(const Vec2i& pos, int fontTex, bool bBack=true)`: Renders the text input in 2D.
-   `onKeyDown(SDL_Event e)`: Handles key press events within the input field.
-   `onText(SDL_Event e)`: Handles text input events within the input field.

### `Commander`

A class that manages a list of commands with key bindings.

-   `commands`: A vector of `Command` objects.
-   `keymap`: A map of key codes to command indices.
-   `add(int id, int key, std::string name="")`: Adds a command to the list.
-   `rebind(int i, int key)`: Rebinds a command to a new key.

### `Command`

A struct that represents a command with an ID, key binding, and name.

-   `id`: The command ID.
-   `key`: The key code for the command.
-   `name`: The command name.

### `GUI_stepper`

A helper class for laying out GUI elements in a vertical column.

-   `x0`, `x1`: The starting and ending y-positions.
-   `step(int n)`: Increments the y-position by a multiple of the default font size.

### `GUIPanelWatcher`

A class that represents a watcher for `GUIPanel` objects.

-   `master`: A pointer to the `GUIPanel` object.
-   `slave`: A pointer to the data value to be watched.
-   `bInt`: A flag indicating whether the value is an integer.
-   `bind(GUIPanel* master_, void* slave_, bool bInt_)`: Binds the watcher to a `GUIPanel` and a data value.
-   `bindLoad(GUIPanel* master_, void* slave_, bool bInt_)`: Binds the watcher and loads the initial value.
-   `apply()`: Applies the value from the `GUIPanel` to the data value.
-   `load()`: Loads the value from the data value to the `GUIPanel`.
-   `check()`: Checks if the `GUIPanel` has been modified and applies the value if necessary.

### `BoundGUI`

A class representing a bound graphical user interface to some set of parameters.

-   `binded`: A flag indicating whether the GUI is bound to data.
-   `drivers`: A pointer to an array of `GUIPanelWatcher` objects.
-   `unbind()`: Unbinds the GUI from the data.
-   `check()`: Checks if any of the bound `GUIPanel` objects have been modified.

## Usage Examples

### Creating a Slider

This example demonstrates how to create a `GUIPanel` as a slider to control a zoom value.

```cpp
#include "GUI.h"

class MyGUI {
public:
    GUI gui;
    float zoom = 10.0f;

    void init() {
        GUI_stepper ylay(1,2);
        // Create a GUIPanel for zoom control
        GUIPanel* zoomPanel = new GUIPanel("Zoom:", 5, ylay.x0, 105, ylay.x1, true, false);
        zoomPanel->setRange(1.0f, 20.0f); // Set the slider range
        zoomPanel->setValue(zoom); // Set the initial value
        zoomPanel->setCommand([this](GUIAbstractPanel* p) { // Set the callback
            zoom = ((GUIPanel*)p)->value; // Update the zoom value
            printf("Zoom: %f\n", zoom);
        });
        gui.addPanel(zoomPanel); // Add the panel to the GUI
    }
};
```

In this example:

-   A `GUI_stepper` is used to manage the vertical layout.
-   A `GUIPanel` is created with the label "Zoom:", positioned at (5, ylay.x0) with a width of 100 and height of ylay.x1.
-   The `isSlider` parameter is set to `true`, making it a slider.
-   The slider's range is set from 1.0 to 20.0.
-   The initial value is set to the `zoom` variable.
-   A lambda function is used as the command callback, which updates the `zoom` variable and prints the new value.

### Creating a Dropdown List

This example demonstrates how to create a `DropDownList` for selecting a view mode.

```cpp
#include "GUI.h"

class MyGUI {
public:
    GUI gui;
    int viewMode = 0;

    void init() {
        GUI_stepper ylay(1,2);
        // Create a DropDownList for view mode selection
        DropDownList* viewModeList = new DropDownList("View Mode:", 5, ylay.x0, 105, 3);
        viewModeList->addItem("Top");
        viewModeList->addItem("Front");
        viewModeList->addItem("Side");
        viewModeList->setCommand([this](GUIAbstractPanel* p) {
            viewMode = ((DropDownList*)p)->iSelected;
            printf("View Mode: %i\n", viewMode);
        });
        gui.addPanel(viewModeList);
    }
};
```

In this example:

-   A `DropDownList` is created with the label "View Mode:", positioned at (5, ylay.x0) with a width of 100 and 3 visible slots.
-   Three items ("Top", "Front", "Side") are added to the list.
-   A lambda function is used as the command callback, which updates the `viewMode` variable and prints the new value.

### Creating a Checkbox List

This example demonstrates how to create a `CheckBoxList` for toggling options.

```cpp
#include "GUI.h"

class MyGUI {
public:
    GUI gui;
    bool showGrid = false;
    bool showAxis = true;

    void init() {
        GUI_stepper ylay(1,2);
        // Create a CheckBoxList for toggling options
        CheckBoxList* optionsList = new CheckBoxList(5, ylay.x0, 105);
        optionsList->caption = "Options";
        optionsList->addBox("Show Grid", &showGrid);
        optionsList->addBox("Show Axis", &showAxis);
        gui.addPanel(optionsList);
    }
};
```

In this example:

-   A `CheckBoxList` is created with the label "Options", positioned at (5, ylay.x0) with a width of 100.
-   Two checkboxes are added: "Show Grid" and "Show Axis", bound to the `showGrid` and `showAxis` variables respectively.
-   The `CheckBoxList` automatically updates the bound variables when the checkboxes are toggled.

These examples provide a basic understanding of how to use the GUI system. You can combine these elements to create more complex and interactive user interfaces.
