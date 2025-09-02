# CriOS Dr. Crowe Coder Platform

## 🚀 Full-Stack Drug Discovery Platform

A modern web platform for the revolutionary CriOS Dr. Crowe Coder system, orchestrating 194 PhD AI agents for accelerated drug discovery.

## Architecture

### Tech Stack

- **Backend**: FastAPI + Python 3.11
- **Frontend**: Next.js 14 + React 18 + TypeScript
- **Database**: PostgreSQL + Redis
- **Real-time**: WebSockets
- **Deployment**: Docker + Kubernetes
- **CI/CD**: GitHub Actions

## Quick Start

### Prerequisites

- Docker & Docker Compose
- Node.js 18+
- Python 3.11+
- PostgreSQL 16+
- Redis 7+

### Local Development

1. **Clone the repository**
```bash
git clone https://github.com/MichaelCrowe11/crios-dr-crowe-coder.git
cd crios-dr-crowe-coder
```

2. **Set up environment variables**
```bash
cp .env.example .env
# Edit .env with your API keys
```

3. **Start with Docker Compose**
```bash
docker-compose up -d
```

4. **Access the platform**
- Frontend: http://localhost:3000
- Backend API: http://localhost:8000
- API Docs: http://localhost:8000/api/docs
- Flower (Celery monitoring): http://localhost:5555

### Manual Setup

#### Backend

```bash
cd platform/backend
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
uvicorn app:app --reload
```

#### Frontend

```bash
cd platform/frontend
npm install
npm run dev
```

## API Endpoints

### Core Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/agents` | List all 194 PhD agents |
| GET | `/api/agents/{id}` | Get specific agent details |
| GET | `/api/pipelines` | List discovery pipelines |
| POST | `/api/discover` | Start discovery pipeline |
| POST | `/api/molecule/analyze` | Analyze molecule |
| POST | `/api/similarity/search` | Similarity search |
| POST | `/api/filter` | Filter compounds |
| POST | `/api/orchestrate` | Custom task orchestration |
| GET | `/api/stats` | System statistics |

### WebSocket Endpoints

- `/ws/pipeline/{id}` - Real-time pipeline updates
- `/ws/agents` - Agent activity stream

## Features

### Dashboard

- **Real-time Monitoring**: Live agent activity and pipeline progress
- **3D Molecule Viewer**: Interactive molecular visualization
- **Agent Network Graph**: D3.js network visualization of 194 agents
- **Pipeline Builder**: Visual pipeline creation and monitoring
- **Statistics**: Performance metrics and success rates

### Discovery Pipelines

1. **Kinase Inhibitors** - Targeted cancer therapies
2. **PROTACs** - Protein degradation
3. **Antibody-Drug Conjugates** - Precision oncology
4. **Natural Products** - Nature-inspired drugs
5. **AI Generative** - ML-designed molecules
6. **Fragment-Based** - Fragment growing/linking
7. **Peptides** - Peptide therapeutics
8. **Drug Repurposing** - New uses for known drugs

## Docker Deployment

### Build Images

```bash
# Build all services
docker-compose build

# Build specific service
docker-compose build backend
docker-compose build frontend
```

### Production Deployment

```bash
# Use production configuration
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d
```

## Testing

### Backend Tests

```bash
cd platform/backend
pytest tests/ -v --cov=src
```

### Frontend Tests

```bash
cd platform/frontend
npm test
npm run test:e2e
```

## CI/CD Pipeline

GitHub Actions workflow automatically:

1. Runs tests on push
2. Builds Docker images
3. Deploys to staging (develop branch)
4. Deploys to production (main branch)
5. Runs security scans
6. Builds documentation

## Environment Variables

### Backend
```env
DATABASE_URL=postgresql://user:pass@localhost/db
REDIS_URL=redis://localhost:6379
ANTHROPIC_API_KEY=your_api_key
SECRET_KEY=your_secret_key
```

### Frontend
```env
NEXT_PUBLIC_API_URL=http://localhost:8000
NEXT_PUBLIC_WS_URL=ws://localhost:8000
```

## Monitoring

- **Prometheus**: http://localhost:9090
- **Grafana**: http://localhost:3001
- **Flower**: http://localhost:5555

## API Documentation

- **Swagger UI**: http://localhost:8000/api/docs
- **ReDoc**: http://localhost:8000/api/redoc

## Contributing

1. Fork the repository
2. Create feature branch (`git checkout -b feature/amazing`)
3. Commit changes (`git commit -m 'Add feature'`)
4. Push to branch (`git push origin feature/amazing`)
5. Open Pull Request

## Security

- Authentication via JWT tokens
- Rate limiting on all endpoints
- Input validation with Pydantic
- SQL injection protection
- XSS protection
- CORS configuration

## Performance

- Redis caching for frequent queries
- Database connection pooling
- Async request handling
- CDN for static assets
- Lazy loading of components

## Scaling

The platform is designed to scale:

- Horizontal scaling with Kubernetes
- Load balancing with Nginx
- Database read replicas
- Redis clustering
- Microservices architecture ready

## License

MIT License - see LICENSE file

## Support

- GitHub Issues: https://github.com/MichaelCrowe11/crios-dr-crowe-coder/issues
- Documentation: https://crios.ai/docs
- Email: support@crios.ai

---

© 2025 Crowe Research. Built with ❤️ for drug discovery.