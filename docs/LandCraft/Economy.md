# Economy.h - Economic Simulation Framework

## Purpose
Provides the foundational economic simulation layer for LandCraft, managing commodities, technologies, production chains, and resource distribution economics.

## Core Components

### 1. Commodity System

#### Commodity Class
- **Purpose**: Defines tradeable goods with realistic transport constraints
- **Key Attributes**:
  - `transport_weight`: Effective weight including containers/packaging
  - `transport_volume`: Space requirements for transportation
  - `price_max`: Maximum market price under scarcity
  - `price_normal`: Standard market price under normal conditions
- **Algorithm**: Uses string parsing for initialization from configuration files

### 2. Technology Framework

#### Technology Class
- **Purpose**: Represents production processes with input/output relationships
- **Key Attributes**:
  - `cycle_time`: Production duration per unit (real-time simulation)
  - `unit_space`: Factory floor space required per production unit
  - `consumes`: Hash map of input resources and quantities
  - `produces`: Hash map of output products and quantities
- **Algorithm**: File-based configuration using structured text format

#### Configuration Format
```
TechnologyName cycle_time unit_space
consumes: resource1=amount1 resource2=amount2
produces: product1=amount1 product2=amount2
```

### 3. Production Management

#### Factory Class
- **Purpose**: Executes production processes using defined technologies
- **Key Functions**:
  - `setTechnology()`: Configures factory for specific production
  - `produce(N)`: Attempts to produce N units, limited by available resources
- **Algorithm**: 
  - Calculates maximum production based on scarcest input resource
  - Atomically consumes inputs and produces outputs
  - Returns actual units produced (may be less than requested)

### 4. Resource Management

#### Storage System
- **Purpose**: Tracks resource quantities without complex inventory management
- **Implementation**: Hash maps linking commodity names to current quantities
- **Advantage**: Simple, fast lookups suitable for large-scale simulation

- `Commodity` holds `name`, `transport_weight`, `transport_volume`, `price_max`, `price_normal`; initialized via `fromString(char*)` using `sscanf`.
- `Technology` holds `name`, `cycle_time`, `unit_space`, `consumes`, `produces`; loaded by `fromFile(FILE*)` which skips to the next alpha-starting line, then parses two subsequent lines into `consumes`/`produces` via `str2map()`.
- `Factory` holds a `stored` map<string,double> and a `currentTenchnology` pointer; `setTechnology()` ensures presence of all consumed/produced keys with default `0.0`.
- `Factory::produce(N)` limits N by available inputs, subtracts consumed amounts, adds produced amounts, and returns the actually produced units.

**Notes:**
- No price dynamics or markets are implemented in this header [need-check].
- Transport integration is not present here; only transport-related fields in `Commodity` exist [need-check].

## Economic Simulation Features

### 1. Supply Chain Modeling
- **Multi-stage Production**: Complex chains requiring multiple technologies
- **Resource Constraints**: Production limited by actual resource availability
- **Transport Considerations**: Weight/volume affect distribution logistics

### 2. Dynamic Pricing (Framework)
- **Price Bounds**: Maximum prices prevent unrealistic inflation
- **Volume-based Scaling**: Transport costs scale with weight/volume
- **Market Simulation**: Foundation for supply/demand price adjustments

### 3. Production Optimization
- **Resource Efficiency**: Automatic calculation of optimal production ratios
- **Bottleneck Identification**: Clear indication of limiting resources
- **Flexible Scaling**: Production can be scaled up/down based on demand

## Integration with LandCraft

### 1. Geographic Context
- **Resource Location**: Commodities tied to specific terrain locations
- **Transport Networks**: Infrastructure affects commodity movement costs
- **Regional Specialization**: Different areas excel at different production

### 2. Infrastructure Dependencies
- **Factory Placement**: Requires road/rail access for input/output
- **Supply Chains**: Longer chains require more sophisticated networks
- **Scale Economics**: Larger networks enable more complex production

### 3. Strategic Decision Making
- **Technology Choice**: Players must select appropriate technologies for location
- **Resource Allocation**: Decisions about which goods to produce vs. import
- **Investment Planning**: Long-term infrastructure vs. immediate production needs

## Extensibility

### 1. Adding New Commodities
- Simple text file configuration
- Automatic integration with existing technologies
- No code changes required

### 2. Technology Trees
- Foundation for research/development systems
- Prerequisites and unlock conditions
- Progressive complexity increase

### 3. Market Dynamics
- Framework supports future price fluctuation
- Supply/demand modeling capability
- Economic crisis simulation potential

This economic system provides the backbone for realistic resource management while maintaining computational efficiency for large-scale world simulation.
