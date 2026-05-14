import json 
import math

class Node:
    def __init__(self, id:str, scatter_coeff:float, sensor:float, threshold: FloatingPointError) -> None:
        self.id = id
        self.scatter_coeff = scatter_coeff
        self.sensor = sensor
        self.threshold = threshold
    
    def __repr__(self) -> str:
        return self.id

class Edge:
    def __init__(self, u: str, v: str, length_m: int, material: str, attenuation: float, noise_coeff, directed: int) -> None:
        self.u = u
        self.v = v
        self.length_m = length_m
        self.material = material
        self.attenuation = attenuation
        self.noise_coeff = noise_coeff
        self.directed = directed
    
    def __repr__(self) -> str:
        return f"{self.u} - {self.v}"
    
    def travel_time_ms(self, material_speeds: dict[str, float]) -> float:
        material = material_speeds[self.material]
        t = self.length_m / material * 1000
        return t

    def other(self, node: Node) -> Node:
        if node.id != self.u and node.id != self.v:
            return None
        if node.id == self.u:
            return self.v
        if not self.directed:
            if node.id == self.v:
                return self.u
        else:
            if node.id == self.v:
                return None
    
    def propagate_energy(self, E_in: float) -> float:
        E_out = E_in * math.exp(-self.attenuation * self.length_m)
        return E_out
    
    def propagate_noise(self, noise_in: float) -> float:
        noise_out = noise_in + self.noise_coeff 
        return noise_out

class SignalSystem:
    def __init__(self) -> None:
        self.vertices = None
        self.edges = None
        self.dct = None
        self.materials = None

    def load_from_json(self, path: str) -> None:
        with open(path, "r", encoding = "utf-8") as file:
            data = json.load(file)
            vertices = data["vertices"]
            edgs = data["edges"]
            self.materials = data["materials"]
            nodes = {}
            edges = {}
            dct = {}
            for item in vertices:
                node = Node(item["id"], item["scatter_coeff"], item["sensor"], item["threshold"])
                nodes[node.id] = node
                dct[node] = []
            for item in edgs:
                edge = Edge(item["u"], item["v"], item["length_m"], item["material"], item["attenuation"], item["noise_coeff"], item["directed"])
                edges[(item["u"], item["v"])] = edge
            
            for ver in nodes.values():
                for edg in edges.values():
                    if ver.id == edg.u:
                        v = nodes[edg.v]
                        dct[ver].append(v)
                        if not edg.directed:
                            dct[v].append(ver)

        self.vertices = nodes
        self.edges = edges
        self.dct = dct

    def dijkstra(self, start: str, R_ms: float, initial_energy: float, alpha: float, beta: float, gamma: float) -> None:
        unvisited = set()
        costs = {}
        times = {}
        energies = {}
        noises = {}
        new_dct = {}
        parents = {}
        for item in self.vertices.values():
            unvisited.add(item)
            costs[item] = float("inf")
            if item.id == start:
                costs[item] = 0
                energies[item] = initial_energy
                times[item] = 0
                noises[item] = 0
        
        while unvisited:
            current = None
            min_cost = float("inf")
            for ver in unvisited:
                if costs[ver] < min_cost:
                    min_cost = costs[ver]
                    current = ver
            if current is not None:
                unvisited.remove(current)
                nbds = self.dct[current]
                for ver in nbds:
                    if ver not in unvisited:
                        continue
                    edge = self.edges.get((current.id, ver.id))
                    if not edge:
                        edge = self.edges.get((ver.id, current.id))
                    time_ms = edge.travel_time_ms(self.materials)
                    time = times[current] + time_ms
                    energy_after_edge = edge.propagate_energy(energies[current])
                    energy = energy_after_edge * (1 - ver.scatter_coeff)
                    noise = edge.propagate_noise(noises[current])
                    cost = gamma * time + alpha*noise + beta * (1/energy)
                    if energy <= 0 or time >= R_ms or cost < 0:
                        continue
                    if cost < costs[ver]:
                        costs[ver] = cost
                        times[ver] = time
                        energies[ver] = energy
                        noises[ver] = noise
                        parents[ver] = current 
                new_dct[current] = {"time": round(times[current],2), "energy": round(energies[current],2), "noise": round(noises[current],2), "cost": round(costs[current],2)}
            else:
                break
        with open("itog", "w", encoding= "utf-8") as file:
            file.write(f"R = {R_ms}:")
            file.write("\n")
            for current in new_dct:
                file.write("\n\n")
                file.write(f"{current.id}:")
                file.write("\n")
                dct_current = new_dct[current]
                for item in dct_current:
                    file.write(f"{item} = {dct_current[item]}")
                    file.write("\n")
                if parents.get(current):
                        path = self.restore_path(current, parents)
                        s = ""
                        for i in range(len(path)):
                            if i != len(path) - 1:
                                s += f" {path[i]} -> "
                            else:
                                s += f" {path[i]}"
                                print("\n")
                        file.write(f"parents: {s}")
                        print("\n")    
                

        

    def restore_path(self, ver, parents):
        path = []
        current = ver
        while current is not None:
            path.append(current.id)
            current = parents.get(current)  
        path.reverse()
        return path

           



ss = SignalSystem()
ss.load_from_json("data.json")
with open("data.json", "r", encoding = "utf-8") as file:
    data = json.load(file)
    ss.dijkstra(data["start"], data["R_ms"], data["initial_energy"], data["alpha"], data["beta"], data["gamma"])