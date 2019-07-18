import { HttpClient, HttpHeaders, HttpParams } from '@angular/common/http';
import { Injectable } from '@angular/core';

export interface Resp {
    // to implement
}

@Injectable({
    providedIn: 'root'
})

export class RequestService {

    constructor(private http: HttpClient) {}

    trimmingStep(body: string) {
        const path = 'http://localhost/trimming';
        return this.http.post<Resp>(
            path,
            body,
            {
                headers: new HttpHeaders({
                    'Content-Type': 'application/json',
                })
            });
    }
}
